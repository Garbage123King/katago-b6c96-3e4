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

def load_norm(channel, scale_index, bias_index):
    norm = nn.BatchNorm2d(channel, affine=True)
    scale = torch.tensor(bins[scale_index]["floats"], dtype=torch.float32)
    bias = torch.tensor(bins[bias_index]["floats"], dtype=torch.float32)
    with torch.no_grad():
        norm.weight.copy_(scale)
        norm.bias.copy_(bias)
    norm.eval()
    return norm

def load_conv(index, k, ic, oc):
    w = torch.tensor(bins[index]["floats"], dtype=torch.float32)
    w = w.view(k, k, ic, oc)
    w = w.permute(3, 2, 0, 1).contiguous()
    padding = 0 if k==1 else 1
    conv = nn.Conv2d(ic, oc, k, padding=padding, bias=False) #不要bias
    with torch.no_grad():
        conv.weight.copy_(w)
    conv.eval()
    return conv

def load_linear(index, ic, oc):
    linear = nn.Linear(ic, oc, bias=False)
    w = torch.tensor(bins[index]["floats"], dtype=torch.float32)
    w = w.view(ic,oc).t()
    with torch.no_grad():
        linear.weight.copy_(w)
    return linear

conv0 = load_conv(0, 3, 22, 96)


trunkScratch=torch.tensor(trunkScratch, dtype=torch.float32).unsqueeze(0)
input19=torch.tensor(input19, dtype=torch.float32).unsqueeze(0) # 22x19x19 变 1x22x19x19
conv0.eval()
output0 = conv0(input19)
diff(output0, trunkScratch)


inputGlobal19 = torch.tensor(inputGlobal19, dtype=torch.float32).unsqueeze(0)
linear0 = load_linear(1, 19, 96)
    
inputMatMulOut=torch.tensor(inputMatMulOut, dtype=torch.float32).unsqueeze(0)
linear0.eval()
output1 = linear0(inputGlobal19)
diff(output1, inputMatMulOut)

output2 = output0 + output1.view(-1, 1, 1)

trunkScratch_afterBias = torch.tensor(trunkScratch_afterBias, dtype=torch.float32).unsqueeze(0)
diff(output2, trunkScratch_afterBias)


#block阶段


def ordi(input, index0, index1, index2, index3, index4, index5):
    norm = load_norm(96, index0, index1)

    # 第一个norm
    output = norm(input)
    output = torch.relu(output)

    #第一个conv
    conv = load_conv(index2, 3, 96, 96)
    output = conv(output)

    #第二个norm
    norm = load_norm(96, index3, index4)
    output = norm(output)
    output = torch.relu(output)

    #第二个conv
    conv = load_conv(index5, 3, 96, 96)
    output = conv(output)

    output = input + output

    return output

def gpool(input, index0, index1, index2, index3, index4, index5, index6, index7, index8, index9):
    norm = load_norm(96, index0, index1)
    output = norm(input)
    output = torch.relu(output)
    conv = load_conv(index2, 3, 96, 64)
    output_备用 = conv(output)

    conv = load_conv(index3, 3, 96, 32)
    output = conv(output)
    norm = load_norm(32, index4, index5)
    output = norm(output)
    output = torch.relu(output)
    output = rowsG(output, False, 19)
    linear = load_linear(index6, 96, 64)
    output = linear(output)

    output = output_备用 + output.view(-1, 1, 1)

    norm = load_norm(64, index7, index8)
    output = norm(output)
    output = torch.relu(output)

    conv = load_conv(index9, 3, 64, 96)
    output = conv(output)

    output = input + output

    return output

def rowsG(input_tensor, is_value_head, board_size):
    """
    input_tensor: PyTorch Tensor [32, 19, 19]
    isValueHead: bool, 决定第3个 channel 输出 mean 相关的缩放值还是 max
    board_size: 实际棋盘大小 (通常为 19)
    """
    device = input_tensor.device
    batch_size = input_tensor.shape[1] # 32
    
    # 1. 计算掩码 (Mask)
    # 逻辑：对于 input[i][j]，如果 i, j >= board_size，则视为 padding 区域
    mask = torch.zeros((19, 19), device=device)
    mask[:board_size, :board_size] = 1.0
    
    # 应用掩码到输入
    masked_input = input_tensor * mask
    
    # 2. 计算 Mean (s / div)
    # div 是棋盘格数，例如 19*19=361
    div = float(board_size * board_size)
    s = torch.sum(masked_input, dim=(2, 3)) # 对 H, W 维度求和，得到 [32]
    mean = s / div
    
    # 3. 计算 Max
    # 逻辑：temp = x + (mask - 1.0)。Padding 区域会减去 1.0 确保不被选为 max
    temp = input_tensor + (mask - 1.0)
    max_val, _ = torch.max(temp.view(batch_size, -1), dim=1) # 展平后求每个 batch 的最大值

    # 4. 构建输出 (Output)
    # 预分配输出空间 [32 * 3] -> [96] 或保持维度以供后续使用
    output = torch.zeros(batch_size * 3, device=device)
    
    sqrtdiv = float(board_size)
    scaling_factor = (sqrtdiv - 14.0)
    
    # Channel 0: Mean
    output[0:32] = mean
    
    # Channel 1: Mean * (sqrtdiv - 14) * 0.1
    output[32:64] = mean * scaling_factor * 0.1
    
    # Channel 2: 根据 isValueHead 决定输出
    if is_value_head:
        # value head 逻辑
        val = mean * (scaling_factor * scaling_factor * 0.01 - 0.1)
        output[64:96] = val
    else:
        # policy head 等逻辑输出 max
        output[64:96] = max_val
        
    return output

# block[0] ordi
output = ordi(output2, 3, 4, 5, 8, 9, 10)
afterconv2 = torch.tensor(afterconv2, dtype=torch.float32).unsqueeze(0)
diff(output, afterconv2)

# block[1] ordi
output = ordi(output, 12, 13, 14, 17, 18, 19)

# block[2] gpool
output = gpool(output, 21, 22, 23, 24, 26, 27, 28, 31, 32, 33)
aftergpool = torch.tensor(aftergpool, dtype=torch.float32).unsqueeze(0)
diff(output, aftergpool)

# block[3] ordi
output = ordi(output, 35, 36, 37, 40, 41, 42)

# block[4] gpool
output = gpool(output, 44, 45, 46, 47, 49, 50, 51, 54, 55, 56)

# block[5] ordi
output = ordi(output, 58, 59, 60, 63, 64, 65)

# final
norm = load_norm(96, 67, 68)
output = norm(output)
output_final = torch.relu(output)

#policy

#tree1
conv = load_conv(69, 1, 96, 32)
output1 = conv(output_final) #[96][19][19] -> [32][19][19]

#tree2
conv = load_conv(70, 1, 96, 32)
output2 = conv(output_final)

norm = load_norm(32, 72, 73)
output = norm(output2)
output = torch.relu(output)

outputG = rowsG(output, False, 19)

linear = load_linear(74, 96, 32)
output = linear(outputG) #[96] -> [32]

output = output1 + output.view(-1, 1, 1)

norm = load_norm(32, 76, 77)
output = norm(output)
output = torch.relu(output)

conv = load_conv(78, 1, 32, 2)
output = conv(output)

policy19 = torch.tensor(policy19, dtype=torch.float32).unsqueeze(0)
diff(output, policy19)

linear = load_linear(79, 96, 32)
output = linear(outputG)

adder = torch.tensor(bins[80]["floats"], dtype=torch.float32)
output = output + adder
output = torch.relu(output)

linear = load_linear(81, 32, 2)
output = linear(output)

pass19 = torch.tensor(pass19, dtype=torch.float32).unsqueeze(0)
diff(output, pass19)

conv = load_conv(82, 1, 96, 32)
output = conv(output_final)
norm = load_norm(32, 84, 85)
output = norm(output)
output_prepare = torch.relu(output)

outputG = rowsG(output_prepare, True, 19)
linear = load_linear(86, 96, 64)
output = linear(outputG)

adder = torch.tensor(bins[87]["floats"], dtype=torch.float32)
output = output + adder
output_bak = torch.relu(output)

linear = load_linear(88, 64, 3)
output = linear(output_bak)

adder = torch.tensor(bins[89]["floats"], dtype=torch.float32)
output = output + adder

value19 = torch.tensor(value19, dtype=torch.float32).unsqueeze(0)
diff(output, value19)

linear = load_linear(90, 64, 6)
output = linear(output_bak)

adder = torch.tensor(bins[91]["floats"], dtype=torch.float32)
output = output + adder

conv = load_conv(92, 1, 32, 1)
output = conv(output_prepare)

ownership = torch.tensor(ownership, dtype=torch.float32).unsqueeze(0)
diff(output, ownership)

'''
    float kernel15[32][96];

    n=0;

    for(int j=0; j <96; j++)
      for(int i=0; i <32; i++)
      {
        kernel15[i][j] = BINS[69].floats[n++];
      }

    float output10[32][19][19];
    conv1x1((float*)output9, 96, (float*)output10, 32, 19, (float*)kernel15);

    /*分支2*/

    float kernel16[32][96];

    n=0;

    for(int j=0; j <96; j++)
      for(int i=0; i <32; i++)
      {
        kernel16[i][j] = BINS[70].floats[n++];
      }

    float output11[32][19][19];
    conv1x1((float*)output9, 96, (float*)output11, 32, 19, (float*)kernel16);

    float scale15[32];
    float bias15[32];
    float output12[32][19][19];
    n=0;

    for(int i=0; i < 32; i++)
    {
      scale15[i] = BINS[72].floats[n];
      bias15[i] = BINS[73].floats[n];
      n++;
    }

    norm(output11, output12, 32, scale15, bias15, board_size);

    float output13[96];
    rowsG(output12, output13, FALSE, board_size);

    // 从这起开始分两个头，policy和pass

    float nn2[32][96];
    n=0;

    for(int j=0; j < 96; j++)
      for(int i=0; i < 32; i++)
      {
        nn2[i][j] = BINS[74].floats[n];
        n++;
      }

    float output14[32];
    conv1x1((float*)output13, 96, (float*)output14, 32, 1, (float*)nn2);

    float output15[32][19][19];
    add_broadcast(output10, output15, output14, 32);

    // 汇聚后再来一次normconv

    float scale16[32];
    float bias16[32];
    n=0;

    for(int i=0; i < 32; i++)
    {
      scale16[i] = BINS[76].floats[n];
      bias16[i] = BINS[77].floats[n];
      n++;
    }
    float output16[32][19][19];

    norm(output15, output16, 32, scale16, bias16, board_size);

    float kernel17[2][32];

    n=0;

    for(int j=0; j <32; j++)
      for(int i=0; i <2; i++)
      {
        kernel17[i][j] = BINS[78].floats[n++];
      }

    conv1x1((float*)output16, 32, (float*)policy, 2, 19, (float*)kernel17);
'''