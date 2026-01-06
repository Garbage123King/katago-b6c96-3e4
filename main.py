import torch
import gzip
import shutil
import struct

gzfile = "model3e4.bin.gz"
binfile = "model3e4.bin"

with gzip.open(gzfile, "rb") as fin, open(binfile, "wb") as fout:
    shutil.copyfileobj(fin, fout)

print("unzip finished:", binfile)

float_len = [19008, 1824, 96, 96, 96, 82944, 96, 96, 96, 96, 82944, 96, 96, 96, 82944, 96, 96, 96, 96, 
82944, 96, 96, 96, 55296, 27648, 32, 32, 32, 6144, 64, 64, 64, 64, 55296, 96, 96, 96, 82944, 96, 96, 96, 
96, 82944, 96, 96, 96, 55296, 27648, 32, 32, 32, 6144, 64, 64, 64, 64, 55296, 96, 96, 96, 82944, 96, 96, 
96, 96, 82944, 96, 96, 96, 3072, 3072, 32, 32, 32, 3072, 32, 32, 32, 64, 3072, 32, 64, 3072, 32, 32, 32, 
6144, 64, 192, 3, 384, 6, 32
]

total_floats = 0
total_min = 1e9
total_max = -1e9

bins=[]

def bytes_to_floats(buf):
    assert len(buf) % 4 == 0
    n = len(buf) // 4
    return list(struct.unpack(f"<{n}f", buf))

with open(binfile, "rb") as fp:
    fp.seek(0)
    raw = fp.read(100000000)
    cnt = raw.count(b"@BIN@")
    raw_split = raw.split(b"@BIN@")
    raw_split_bin = raw_split[1:]
    raw_split_description = raw_split[0:-1]
    assert(cnt == len(float_len))
    assert(cnt == len(raw_split_bin))
    assert(cnt == len(raw_split_description))
    for i in range(cnt):
        bin = {}
        if i==0:
            bin["description"] = raw_split_description[i].decode("ascii").replace('\n',' ')
        else:
            bin["description"] = raw_split_description[i][float_len[i-1]*4:].decode("ascii").replace('\n',' ')

        bin["binary"] = raw_split_bin[i][0:float_len[i]*4]
        print(bin["description"], ", binary lenth: ", len(bin["binary"]), ", theory: ", float_len[i]*4)
        bin["floats"]=bytes_to_floats(bin["binary"])
        bins.append(bin)


