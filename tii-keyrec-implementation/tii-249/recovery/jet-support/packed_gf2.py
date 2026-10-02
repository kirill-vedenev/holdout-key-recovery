"""Use Sage's native bit-packed PNG deserializer for integer bit rows.

Sage 10.7 serializes GF(2) matrices as indexed PNGs. Constructing that format
avoids the Python-level conversion of each individual matrix entry.
"""
import struct
import zlib
from sage.all import GF, matrix
from sage.matrix.matrix_mod2_dense import unpickle_matrix_mod2_dense_v2

_translate = bytes(255-int(f'{i:08b}'[::-1], 2) for i in range(256))


def _chunk(kind, data):
    return struct.pack('>I', len(data))+kind+data+struct.pack('>I', zlib.crc32(kind+data))


def packed_matrix(rows, width, max_chunk_entries=64*1024*1024):
    count = len(rows)
    if not count or not width:
        return matrix(GF(2), count, width)
    # libgd temporarily allocates one byte per pixel and rejects products
    # exceeding INT_MAX. Bounded chunks also keep its transient RAM small.
    if count*width > max_chunk_entries and count > 1:
        height = max(1, max_chunk_entries//width)
        blocks = [packed_matrix(rows[a:a+height],width,max_chunk_entries)
                  for a in range(0,count,height)]
        while len(blocks) > 1:
            blocks = [blocks[a].stack(blocks[a+1]) if a+1 < len(blocks) else blocks[a]
                      for a in range(0,len(blocks),2)]
        return blocks[0]
    stride = (width+7)//8
    scanlines = []
    for value in rows:
        value = int(value)
        if value < 0 or value.bit_length() > width:
            raise ValueError('Input row exceeds the declared width')
        # PNG stores bits most significant first, and Sage encodes 1 as
        # black palette index 0, hence bit reversal and complementation.
        scanlines.append(b'\0'+value.to_bytes(stride, 'little').translate(_translate))
    png = (b'\x89PNG\r\n\x1a\n'
           + _chunk(b'IHDR', struct.pack('>IIBBBBB', width, count, 1, 3, 0, 0, 0))
           + _chunk(b'PLTE', b'\0\0\0\xff\xff\xff')
           + _chunk(b'IDAT', zlib.compress(b''.join(scanlines), 1))
           + _chunk(b'IEND', b''))
    return unpickle_matrix_mod2_dense_v2(count, width, png, len(png), False)
