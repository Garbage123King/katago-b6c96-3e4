import torch
import gzip
import shutil
import struct
import torch.nn as nn
from checkresult import *

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

def diff(output0, output1):
    diffmax = torch.max(torch.abs(output0 - output1))
    diffsum = torch.sum(torch.abs(output0 - output1))
    print("diffmax: ", diffmax.item(), "diffsum: ", diffsum.item())

device = "cpu"

w = torch.tensor(bins[0]["floats"], dtype=torch.float32)
w = w.view(3,3,22,96)
w = w.permute(3, 2, 0, 1).contiguous()
conv0 = nn.Conv2d(22, 96, 3, padding=1, bias=False) #不要bias
with torch.no_grad():
    conv0.weight.copy_(w)


trunkScratch=torch.tensor(trunkScratch, dtype=torch.float32).unsqueeze(0)
input19=torch.tensor(input19, dtype=torch.float32).unsqueeze(0) # 22x19x19 变 1x22x19x19
conv0.eval()
output0 = conv0(input19)
diff(output0, trunkScratch)


inputGlobal19 = torch.tensor(inputGlobal19, dtype=torch.float32).unsqueeze(0)
linear0 = nn.Linear(19, 96, bias=False)
w = torch.tensor(bins[1]["floats"], dtype=torch.float32)
w = w.view(19,96).t()
with torch.no_grad():
    linear0.weight.copy_(w)
    
inputMatMulOut=torch.tensor(inputMatMulOut, dtype=torch.float32).unsqueeze(0)
linear0.eval()
output1 = linear0(inputGlobal19)
diff(output1, inputMatMulOut)

output2 = output0 + output1.view(-1, 1, 1)

trunkScratch_afterBias = torch.tensor(trunkScratch_afterBias, dtype=torch.float32).unsqueeze(0)
diff(output2, trunkScratch_afterBias)

#block阶段

norm0 = nn.BatchNorm2d(96, affine=True)
scale0 = torch.tensor(bins[3]["floats"], dtype=torch.float32)
bias0 = torch.tensor(bins[4]["floats"], dtype=torch.float32)
with torch.no_grad():
    norm0.weight.copy_(scale0)
    norm0.bias.copy_(bias0)

afternorm=torch.tensor(afternorm, dtype=torch.float32).unsqueeze(0)
norm0.eval()
output3 = norm0(output2)
output3 = torch.relu(output3)
diff(output3, afternorm)


'''
norm(input, output0, 96, scale0, bias0, board_size);

conv3x3(output0, 96, output1, 96, (float (*)[3][3])kernel1);

norm(output1, output2, 96, scale1, bias1, board_size);

conv3x3(output2, 96, output3, 96, (float (*)[3][3])kernel2);

add(input, output, output3);
'''