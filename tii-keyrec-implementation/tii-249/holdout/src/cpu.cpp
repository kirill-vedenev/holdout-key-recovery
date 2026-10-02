#include "engine.hpp"
namespace tii {
struct CpuBackend final:Backend {
    std::vector<u64> a;
    explicit CpuBackend(const Request& request):Backend(request),a(r.rows.size()*r.words) {}
    std::string name() const override { return "cpu_exact_panel_elimination"; }
    void initialize() override { a=assemble_cpu(r); }
    void gather(u64 word,std::vector<u64>& columns) override {
        columns.resize(r.rows.size());for(u64 i=0;i<r.rows.size();++i) columns[i]=a[i*r.words+word];
    }
    void read_rows(u64 first,u64 count,u64* dst) override { std::copy_n(a.data()+first*r.words,count*r.words,dst); }
    void write_rows(u64 first,u64 count,const u64* src) override { std::copy_n(src,count*r.words,a.data()+first*r.words); }
    void eliminate(u64 begin,const std::vector<u32>& pivots,const std::vector<u64>& panel) override {
        const u64 start=begin/64,tail=r.words-start,entries=u64(1)<<pivots.size();
        std::vector<u64> table(entries*tail);
        for(u64 j=1;j<entries;++j) {
            int bit=__builtin_ctzll(j);u64 prev=j&(j-1);
            for(u64 w=0;w<tail;++w) table[j*tail+w]=table[prev*tail+w]^panel[u64(bit)*r.words+start+w];
        }
        for(u64 row=0;row<r.rows.size();++row) {
            unsigned index=0;
            for(unsigned j=0;j<pivots.size();++j) index|=((a[row*r.words+pivots[j]/64]>>(pivots[j]%64))&1)<<j;
            if(index) for(u64 w=0;w<tail;++w) a[row*r.words+start+w]^=table[u64(index)*tail+w];
        }
    }
    std::vector<u64> extract(const std::vector<u64>& values,const std::vector<u64>& free_mask,
                             const std::vector<u32>& row_pivot,u32 count) override {
        u64 sw=(u64(count)+63)/64;auto result=values;
        for(u64 row=0;row<r.rows.size();++row) if(row_pivot[row]!=NONE) {
            u64 dest=u64(row_pivot[row])*sw;
            for(u64 w=0;w<r.words;++w) {
                u64 bits=a[row*r.words+w]&free_mask[w];
                while(bits) {
                    u64 c=w*64+__builtin_ctzll(bits);bits&=bits-1;
                    for(u64 j=0;j<sw;++j) result[dest+j]^=values[c*sw+j];
                }
            }
        }
        return result;
    }
};
}
int main(int argc,char** argv) {
    return tii::cli(argc,argv,[](const tii::Request& r,const tii::Options&) {
        return std::make_unique<tii::CpuBackend>(r);
    });
}
