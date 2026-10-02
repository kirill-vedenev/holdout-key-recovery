#pragma once
#include "common.hpp"
#include <memory>

namespace tii {
struct Options {
    u32 count=0;
    int panel_bits=8,threads=16;
    u64 seed=17320260930ULL,checkpoint_every=2048,max_panels=0;
    double max_gib=64;
    bool resume=false;
    std::vector<int> devices{0,1};
};
struct Backend {
    const Request& r;
    explicit Backend(const Request& request):r(request) {}
    virtual ~Backend()=default;
    virtual std::string name() const=0;
    virtual void initialize()=0;
    virtual void gather(u64 word,std::vector<u64>& columns)=0;
    virtual void read_rows(u64 first,u64 count,u64* dst)=0;
    virtual void write_rows(u64 first,u64 count,const u64* src)=0;
    virtual void eliminate(u64 begin,const std::vector<u32>& pivots,const std::vector<u64>& panel)=0;
    virtual std::vector<u64> extract(const std::vector<u64>& free_values,const std::vector<u64>& free_mask,
                                    const std::vector<u32>& row_pivot,u32 count)=0;
};
inline void save_checkpoint(Backend& b,const fs::path& path,u64 next,const std::vector<u32>& pivots) {
    fs::path tmp=path.string()+".tmp"; std::ofstream f(tmp,std::ios::binary);
    require(bool(f),"cannot create checkpoint"); f.write("TCP1",4); f.write(b.r.id.data(),64);
    put64(f,b.r.rows.size());put64(f,b.r.N);put64(f,b.r.words);put64(f,next);
    for(u32 p:pivots) put32(f,p);
    const u64 chunk=std::max<u64>(1,(u64(64)<<20)/(b.r.words*8));
    std::vector<u64> buffer(std::min<u64>(chunk,b.r.rows.size())*b.r.words);
    for(u64 first=0;first<b.r.rows.size();first+=chunk) {
        u64 count=std::min<u64>(chunk,b.r.rows.size()-first);
        b.read_rows(first,count,buffer.data());
        f.write(reinterpret_cast<char*>(buffer.data()),count*b.r.words*8);
    }
    finish_file(f,tmp,path);
    std::cout<<"{\"event\":\"checkpoint\",\"next_column\":"<<next<<"}"<<std::endl;
}
inline u64 load_checkpoint(Backend& b,const fs::path& path,std::vector<u32>& pivots) {
    std::ifstream f(path,std::ios::binary); require(bool(f),"resume checkpoint not found");
    char magic[4],id[64]; f.read(magic,4);f.read(id,64);
    require(f && std::memcmp(magic,"TCP1",4)==0 && std::string(id,64)==b.r.id,"checkpoint belongs to another instance");
    u64 rows=get_number(f,8),N=get_number(f,8),words=get_number(f,8),next=get_number(f,8);
    require(rows==b.r.rows.size() && N==b.r.N && words==b.r.words && next<=N,"checkpoint dimensions invalid");
    require(fs::file_size(path)==100+rows*4+rows*words*8,"checkpoint file length mismatch");
    std::vector<bool> seen(N,false);
    for(u64 i=0;i<rows;++i) {
        pivots[i]=get_number(f,4);
        if(pivots[i]!=NONE) {
            require(pivots[i]<next && !seen[pivots[i]],"checkpoint pivot metadata invalid");
            seen[pivots[i]]=true;
        }
    }
    const u64 chunk=std::max<u64>(1,(u64(64)<<20)/(words*8));
    std::vector<u64> buffer(std::min<u64>(chunk,rows)*words);
    for(u64 first=0;first<rows;first+=chunk) {
        u64 count=std::min(chunk,rows-first);
        f.read(reinterpret_cast<char*>(buffer.data()),count*words*8);require(bool(f),"checkpoint matrix truncated");
        b.write_rows(first,count,buffer.data());
    }
    return next;
}
inline u64 random64(u64& state) {
    u64 z=(state+=0x9e3779b97f4a7c15ULL);
    z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;
    z=(z^(z>>27))*0x94d049bb133111ebULL;
    return z^(z>>31);
}
inline int solve(Backend& b,const fs::path& out,const Options& options) {
    const Request& r=b.r; Options opt=options;
    if(!opt.count) opt.count=(r.minimum+63)/64*64;
    require(opt.count>=u32(r.minimum) && opt.count<=r.N,"sample count violates the requested minimum or ambient dimension");
    require(opt.panel_bits>=1 && opt.panel_bits<=8,"panel bits must be in [1,8]");
    require(!fs::exists(out/"solver.json") && !fs::exists(out/"kernel-colex.bin"),"solver output already exists");
    fs::create_directories(out);
    std::signal(SIGINT,on_signal);std::signal(SIGTERM,on_signal);
    double start=now(); std::vector<u32> row_pivot(r.rows.size(),NONE);
    u64 column=0;
    if(opt.resume) column=load_checkpoint(b,out/"matrix.checkpoint",row_pivot);
    else {
        require(!fs::exists(out/"matrix.checkpoint"),"checkpoint exists; use --resume");
        b.initialize();
    }
    double assembled=now();
    std::vector<u32> active; std::vector<bool> pivot_column(r.N,false);u64 rank=0;
    for(u32 row=0;row<row_pivot.size();++row) {
        if(row_pivot[row]==NONE) active.push_back(row);
        else { pivot_column[row_pivot[row]]=true; ++rank; }
    }
    std::cout<<"{\"event\":\"elimination_start\",\"backend\":\""<<escape(b.name())
             <<"\",\"rows\":"<<r.rows.size()<<",\"columns\":"<<r.N
             <<",\"matrix_bytes\":"<<r.rows.size()*r.words*8<<",\"next_column\":"<<column<<"}"<<std::endl;
    std::vector<u64> gathered; u64 panels=0; double last=now();
    while(column<r.N && !active.empty()) {
        if(interrupted || (opt.max_panels && panels>=opt.max_panels)) {
            save_checkpoint(b,out/"matrix.checkpoint",column,row_pivot);
            std::cout<<"{\"event\":\"paused\",\"reason\":\""<<(interrupted?"signal":"panel_limit")<<"\"}"<<std::endl;
            return 75;
        }
        const int width=std::min<u64>({u64(opt.panel_bits),64-column%64,r.N-column});
        b.gather(column/64,gathered);
        u64 local_basis[8]{};std::vector<u32> chosen;
        std::size_t i=0;
        while(i<active.size() && chosen.size()<std::size_t(width)) {
            u32 row=active[i]; u64 mask=(gathered[row]>>(column%64))&((u64(1)<<width)-1);
            while(mask) {
                int bit=__builtin_ctzll(mask);
                if(local_basis[bit]) mask^=local_basis[bit];
                else { local_basis[bit]=mask;break; }
            }
            if(mask) { chosen.push_back(row); active[i]=active.back();active.pop_back(); }
            else ++i;
        }
        if(!chosen.empty()) {
            std::vector<u64> panel(chosen.size()*r.words);
            for(std::size_t a=0;a<chosen.size();++a) b.read_rows(chosen[a],1,panel.data()+a*r.words);
            std::vector<u32> pivots;std::size_t lead=0;
            for(int bit=0;bit<width && lead<chosen.size();++bit) {
                u64 col=column+bit;std::size_t found=lead;
                while(found<chosen.size() && !((panel[found*r.words+col/64]>>(col%64))&1)) ++found;
                if(found==chosen.size()) continue;
                if(found!=lead) for(u64 w=0;w<r.words;++w) std::swap(panel[lead*r.words+w],panel[found*r.words+w]);
                for(std::size_t a=0;a<chosen.size();++a) if(a!=lead && ((panel[a*r.words+col/64]>>(col%64))&1))
                    for(u64 w=column/64;w<r.words;++w) panel[a*r.words+w]^=panel[lead*r.words+w];
                pivots.push_back(col);++lead;
            }
            require(lead==chosen.size(),"panel rank inconsistency");
            b.eliminate(column,pivots,panel);
            for(std::size_t a=0;a<chosen.size();++a) {
                b.write_rows(chosen[a],1,panel.data()+a*r.words);
                row_pivot[chosen[a]]=pivots[a];pivot_column[pivots[a]]=true;
            }
            rank+=chosen.size();
        }
        column+=width;++panels;
        if(now()-last>=5 || column==r.N) {
            std::cout<<"{\"event\":\"progress\",\"columns_processed\":"<<column
                     <<",\"rank_so_far\":"<<rank<<",\"elapsed_seconds\":"<<now()-start<<"}"<<std::endl;last=now();
        }
        if(opt.checkpoint_every && panels%opt.checkpoint_every==0)
            save_checkpoint(b,out/"matrix.checkpoint",column,row_pivot);
    }
    column=r.N;double eliminated=now();
    if(interrupted) { save_checkpoint(b,out/"matrix.checkpoint",column,row_pivot); return 75; }
    std::vector<u32> free;
    for(u32 c=0;c<r.N;++c) if(!pivot_column[c]) free.push_back(c);
    require(free.size()>=opt.count,"actual kernel is too small for the requested independent sample");
    const u64 sw=(u64(opt.count)+63)/64;
    std::vector<u64> free_mask(r.words),values(r.N*sw);u64 state=opt.seed;
    for(u32 c:free) {
        free_mask[c/64]|=u64(1)<<(c%64);
        for(u64 w=0;w<sw;++w) values[u64(c)*sw+w]=random64(state);
        if(opt.count%64) values[u64(c)*sw+sw-1]&=(u64(1)<<(opt.count%64))-1;
    }
    // A random choice of free coordinates carries an identity minor. Other free coordinates
    // are dense seeded values. This guarantees independence without claiming uniform sampling.
    for(std::size_t i=free.size();i>1;--i) std::swap(free[i-1],free[random64(state)%i]);
    for(u32 a=0;a<opt.count;++a) {
        u64 c=free[a];std::fill_n(values.data()+c*sw,sw,0);values[c*sw+a/64]=u64(1)<<(a%64);
    }
    std::cout<<"{\"event\":\"extract\",\"rank\":"<<rank<<",\"nullity\":"<<free.size()
             <<",\"elements\":"<<opt.count<<"}"<<std::endl;
    auto columns=b.extract(values,free_mask,row_pivot,opt.count);
    write_kernel(out/"kernel-colex.bin",r,opt.count,columns);
    std::ostringstream report; report<<std::setprecision(12)
        <<"{\"schema\":\"tii-holdout-solver-v1\",\"instance_id\":\""<<r.id
        <<"\",\"backend\":\""<<escape(b.name())<<"\",\"heldout_original\":"<<r.heldout_original
        <<",\"single_position\":true,\"observed_rank\":"<<rank<<",\"observed_nullity\":"<<free.size()
        <<",\"elements\":"<<opt.count<<",\"minimum_elements\":"<<r.minimum
        <<",\"seed\":"<<opt.seed<<",\"sampling\":\"random_free_values_with_identity_minor\""
        <<",\"initialization_seconds\":"<<assembled-start<<",\"elimination_seconds_this_invocation\":"<<eliminated-assembled
        <<",\"extraction_seconds\":"<<now()-eliminated<<",\"membership_verified\":false,\"identity_free_coordinates\":[";
    for(u32 a=0;a<opt.count;++a) { if(a) report<<',';report<<free[a]; } report<<"]}";
    json_file(out/"solver.json",report.str());std::cout<<report.str()<<std::endl;
    return 0;
}
inline u64 unsigned_arg(const std::string& s) {
    require(!s.empty() && s.find_first_not_of("0123456789")==std::string::npos,"expected unsigned integer");
    return std::stoull(s);
}
inline Options parse_options(int argc,char** argv,int start) {
    Options o;
    for(int i=start;i<argc;++i) {
        std::string k=argv[i];
        if(k=="--resume") { o.resume=true;continue; }
        require(i+1<argc,"missing argument for "+k);std::string v=argv[++i];
        if(k=="--count") { u64 x=unsigned_arg(v);require(x<=NONE,"count too large");o.count=x; }
        else if(k=="--seed") o.seed=unsigned_arg(v);
        else if(k=="--panel-bits") {u64 x=unsigned_arg(v);require(x>=1 && x<=8,"panel bits outside [1,8]");o.panel_bits=x;}
        else if(k=="--threads") {u64 x=unsigned_arg(v);require(x>=1 && x<=256,"threads outside [1,256]");o.threads=x;}
        else if(k=="--checkpoint-every") o.checkpoint_every=unsigned_arg(v);
        else if(k=="--max-panels") o.max_panels=unsigned_arg(v);
        else if(k=="--max-gib") { std::size_t p; o.max_gib=std::stod(v,&p);require(p==v.size() && o.max_gib>0 && o.max_gib<=1024,"invalid memory cap"); }
        else if(k=="--devices") {
            o.devices.clear();std::istringstream in(v);std::string d;
            while(std::getline(in,d,',')) {u64 x=unsigned_arg(d);require(x<1024,"device ID out of range");o.devices.push_back(x);}
            require(!o.devices.empty() && o.devices.size()<=16,"invalid device count");
            auto copy=o.devices;std::sort(copy.begin(),copy.end());require(std::adjacent_find(copy.begin(),copy.end())==copy.end(),"duplicate GPU IDs");
        } else throw std::runtime_error("unknown option: "+k);
    }
    return o;
}
template<class Factory> int cli(int argc,char** argv,Factory factory) {
    try {
        u32 endian=1;require(*reinterpret_cast<unsigned char*>(&endian)==1,"only little-endian hosts are supported");
        if(argc>=3 && std::string(argv[1])=="inspect") {
            auto r=read_request(argv[2]);
            std::cout<<"{\"rows\":"<<r.rows.size()<<",\"columns\":"<<r.N<<",\"matrix_bytes\":"<<r.rows.size()*r.words*8
                     <<",\"minimum_elements\":"<<r.minimum<<",\"heldout_original\":"<<r.heldout_original<<"}"<<std::endl;return 0;
        }
        if(argc>=5 && std::string(argv[1])=="verify") {
            auto r=read_request(argv[2]);auto k=read_kernel(argv[3]);auto o=parse_options(argc,argv,5);
            require(!fs::exists(argv[4]),"verification output already exists");auto report=verify(r,k,o.threads);
            json_file(argv[4],report);std::cout<<report<<std::endl;return 0;
        }
        if(argc==4 && std::string(argv[1])=="export") {export_whk(read_kernel(argv[2]),argv[3]);return 0;}
        if(argc>=4 && std::string(argv[1])=="solve") {
            auto r=read_request(argv[2]);auto o=parse_options(argc,argv,4);
            require(double(r.rows.size())*r.words*8<=o.max_gib*1073741824.,"matrix exceeds configured total memory cap");
            auto backend=factory(r,o);return solve(*backend,argv[3],o);
        }
        std::cerr<<"usage: holdout-{cpu,cuda} inspect OPERATOR | solve OPERATOR OUTDIR [--count N --devices 0,1 --resume --max-panels N --checkpoint-every N --panel-bits 8 --max-gib 64 --seed N] | verify OPERATOR KERNEL REPORT [--threads 16] | export KERNEL WHK_OUTPUT\n";
        return 2;
    } catch(const std::exception& e) {std::cerr<<"Experiment refused: "<<e.what()<<std::endl;return 1;}
}
} // namespace tii
