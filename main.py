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
    """全局池化逻辑封装"""
    batch_size, channels, h, w = input_tensor.shape
    device = input_tensor.device
    
    mask = torch.zeros((h, w), device=device)
    mask[:board_size, :board_size] = 1.0
    masked_input = input_tensor * mask
    
    # Mean
    div = float(board_size * board_size)
    mean = torch.sum(masked_input, dim=(2, 3)) / div
    
    # Max
    temp = input_tensor + (mask - 1.0)
    max_val, _ = torch.max(temp.view(batch_size, channels, -1), dim=2)

    sqrtdiv = float(board_size)
    scaling_factor = (sqrtdiv - 14.0)
    
    # 拼接通道: [batch, 32] -> [batch, 96]
    ch1 = mean
    ch2 = mean * scaling_factor * 0.1
    if is_value_head:
        ch3 = mean * (scaling_factor * scaling_factor * 0.01 - 0.1)
    else:
        ch3 = max_val
        
    return torch.cat([ch1, ch2, ch3], dim=1)

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


#class版
print("开始class版")

class OrdiBlock(nn.Module):
    def __init__(self, channels=96):
        super().__init__()
        self.norm1 = nn.BatchNorm2d(channels)
        self.relu1 = nn.ReLU(inplace=True)
        self.conv1 = nn.Conv2d(channels, channels, 3, padding=1, bias=False)
        self.norm2 = nn.BatchNorm2d(channels)
        self.relu2 = nn.ReLU(inplace=True)
        self.conv2 = nn.Conv2d(channels, channels, 3, padding=1, bias=False)
    
        
    def apply_mask(self, tensor, board_size):
        """将棋盘范围外的区域清零"""
        # 假设 tensor 形状为 [batch, channels, height, width]
        # 如果当前输入的 H 或 W 大于 board_size，则执行 mask
        if tensor.shape[2] > board_size or tensor.shape[3] > board_size:
            # 这种切片赋值法在 inplace=True 时比较高效
            tensor[:, :, board_size:, :] = 0
            tensor[:, :, :, board_size:] = 0
        return tensor

    def forward(self, x, board_size):
        identity = x
        out = self.norm1(x)
        out = self.apply_mask(out, board_size) # board_size外的清零
        out = self.relu1(out)
        out = self.conv1(out)
        out = self.norm2(out)
        out = self.apply_mask(out, board_size) # board_size外的清零
        out = self.relu2(out)
        out = self.conv2(out)
        return out + identity

class GPoolBlock(nn.Module):
    def __init__(self, in_ch=96, mid_ch_c=64, mid_ch_g=32):
        super().__init__()
        self.norm1 = nn.BatchNorm2d(in_ch)
        self.relu1 = nn.ReLU(inplace=True)
        
        self.conv_main = nn.Conv2d(in_ch, mid_ch_c, 3, padding=1, bias=False)
        self.conv_gpool = nn.Conv2d(in_ch, mid_ch_g, 3, padding=1, bias=False)
        
        self.norm_g = nn.BatchNorm2d(mid_ch_g)
        self.relu_g = nn.ReLU(inplace=True)
        
        self.linear_g = nn.Linear(mid_ch_g * 3, mid_ch_c, bias=False)
        
        self.norm2 = nn.BatchNorm2d(mid_ch_c)
        self.relu2 = nn.ReLU(inplace=True)
        self.conv_final = nn.Conv2d(mid_ch_c, in_ch, 3, padding=1, bias=False)


    def apply_mask(self, tensor, board_size):
        """将棋盘范围外的区域清零"""
        # 假设 tensor 形状为 [batch, channels, height, width]
        # 如果当前输入的 H 或 W 大于 board_size，则执行 mask
        if tensor.shape[2] > board_size or tensor.shape[3] > board_size:
            # 这种切片赋值法在 inplace=True 时比较高效
            tensor[:, :, board_size:, :] = 0
            tensor[:, :, :, board_size:] = 0
        return tensor

    def forward(self, x, board_size):
        identity = x
        out = self.norm1(x)
        out = self.apply_mask(out, board_size) # board_size外的清零
        out = self.relu1(out)
        
        main_feat = self.conv_main(out)
        
        # GPool 分支
        g = self.conv_gpool(out)
        g = self.norm_g(g)
        g = self.apply_mask(g, board_size) # board_size外的清零
        g = self.relu_g(g)
        g_vec = rowsG(g, False, board_size)
        g_feat = self.linear_g(g_vec)
        
        # 融合
        out = main_feat + g_feat.view(g_feat.size(0), g_feat.size(1), 1, 1)
        out = self.norm2(out)
        out = self.apply_mask(out, board_size) # board_size外的清零
        out = self.relu2(out)
        out = self.conv_final(out)
        return out + identity

class KataNet(nn.Module):
    def __init__(self):
        super().__init__()
        # Input Layer
        self.conv0 = nn.Conv2d(22, 96, 3, padding=1, bias=False)
        self.linear0 = nn.Linear(19, 96, bias=False) # 针对 GlobalInput
        
        # Blocks
        self.layer0 = OrdiBlock(96)
        self.layer1 = OrdiBlock(96)
        self.layer2 = GPoolBlock(96, 64, 32)
        self.layer3 = OrdiBlock(96)
        self.layer4 = GPoolBlock(96, 64, 32)
        self.layer5 = OrdiBlock(96)
        
        # Final Norm
        self.final_norm = nn.BatchNorm2d(96)
        
        # Policy Head
        self.p_conv1 = nn.Conv2d(96, 32, 1, bias=False)
        self.p_conv2 = nn.Conv2d(96, 32, 1, bias=False)
        self.p_norm2 = nn.BatchNorm2d(32)
        self.p_linear_g = nn.Linear(96, 32, bias=False)
        self.p_norm_combine = nn.BatchNorm2d(32)
        self.p_conv_final = nn.Conv2d(32, 2, 1, bias=False)
        
        # Pass & Value & Ownership (此处省略部分重复定义的层，逻辑同上)
        # ... 可以根据需要继续添加 ...

    def apply_mask(self, tensor, board_size):
        """将棋盘范围外的区域清零"""
        # 假设 tensor 形状为 [batch, channels, height, width]
        # 如果当前输入的 H 或 W 大于 board_size，则执行 mask
        if tensor.shape[2] > board_size or tensor.shape[3] > board_size:
            # 这种切片赋值法在 inplace=True 时比较高效
            tensor[:, :, board_size:, :] = 0
            tensor[:, :, :, board_size:] = 0
        return tensor
    
    def forward(self, img_input, global_input, board_size):
        # 初始融合
        x = self.conv0(img_input)
        g_feat = self.linear0(global_input)
        x = x + g_feat.view(g_feat.size(0), g_feat.size(1), 1, 1)
        
        # 特征提取
        x = self.layer0(x, board_size)
        x = self.layer1(x, board_size)
        x = self.layer2(x, board_size)
        x = self.layer3(x, board_size)
        x = self.layer4(x, board_size)
        x = self.layer5(x, board_size)
        
        x = torch.relu(self.final_norm(x))
        
        # Policy 分支 (示例)
        p1 = self.p_conv1(x)
        p2 = torch.relu(self.apply_mask(self.p_norm2(self.p_conv2(x)), board_size))
        p2_g = rowsG(p2, False, board_size)
        p2_feat = self.p_linear_g(p2_g)
        
        policy_out = p1 + p2_feat.view(p2_feat.size(0), p2_feat.size(1), 1, 1)
        policy_out = self.p_norm_combine(policy_out)
        policy_out = self.apply_mask(policy_out, board_size) # board_size外的清零
        policy_out = torch.relu(policy_out)
        final_policy_map = self.p_conv_final(policy_out) 

        return final_policy_map
        
def load_weights_from_bins(model, bins):
    """
    根据你在脚本中的 index 顺序，手动给 model 的各层 weight 赋值。
    例如：
    """
    def set_conv_weight(layer, bin_idx, k, ic, oc):
        w = torch.tensor(bins[bin_idx]["floats"]).view(k, k, ic, oc).permute(3, 2, 0, 1)
        layer.weight.data.copy_(w)

    def set_norm_weight(layer, scale_idx, bias_idx):
        layer.weight.data.copy_(torch.tensor(bins[scale_idx]["floats"]))
        layer.bias.data.copy_(torch.tensor(bins[bias_idx]["floats"]))

    def set_linear_weight(layer, bin_idx, ic, oc):
        w = torch.tensor(bins[bin_idx]["floats"]).view(ic, oc).t()
        layer.weight.data.copy_(w)

    # 按照你的脚本 index 开始赋值
    set_conv_weight(model.conv0, 0, 3, 22, 96)
    set_linear_weight(model.linear0, 1, 19, 96)

    # Layer 0 (Ordi)
    set_norm_weight(model.layer0.norm1, 3, 4)
    set_conv_weight(model.layer0.conv1, 5, 3, 96, 96)
    set_norm_weight(model.layer0.norm2, 8, 9)
    set_conv_weight(model.layer0.conv2, 10, 3, 96, 96)
    
    # Layer 1 (Ordi)
    set_norm_weight(model.layer1.norm1, 12, 13)
    set_conv_weight(model.layer1.conv1, 14, 3, 96, 96)
    set_norm_weight(model.layer1.norm2, 17, 18)
    set_conv_weight(model.layer1.conv2, 19, 3, 96, 96)

    # Layer 2 (Gpool)
    set_norm_weight(model.layer2.norm1, 21, 22)
    set_conv_weight(model.layer2.conv_main, 23, 3, 96, 64)
    set_conv_weight(model.layer2.conv_gpool, 24, 3, 96, 32)
    set_norm_weight(model.layer2.norm_g, 26, 27)
    set_linear_weight(model.layer2.linear_g, 28, 96, 64)
    set_norm_weight(model.layer2.norm2, 31, 32)
    set_conv_weight(model.layer2.conv_final, 33, 3, 64, 96)

    # Layer 3 (Ordi)
    set_norm_weight(model.layer3.norm1, 35, 36)
    set_conv_weight(model.layer3.conv1, 37, 3, 96, 96)
    set_norm_weight(model.layer3.norm2, 40, 41)
    set_conv_weight(model.layer3.conv2, 42, 3, 96, 96)

    # Layer 4 (Gpool)
    set_norm_weight(model.layer4.norm1, 44, 45)
    set_conv_weight(model.layer4.conv_main, 46, 3, 96, 64)
    set_conv_weight(model.layer4.conv_gpool, 47, 3, 96, 32)
    set_norm_weight(model.layer4.norm_g, 49, 50)
    set_linear_weight(model.layer4.linear_g, 51, 96, 64)
    set_norm_weight(model.layer4.norm2, 54, 55)
    set_conv_weight(model.layer4.conv_final, 56, 3, 64, 96)

    # Layer 5 (Ordi)
    set_norm_weight(model.layer5.norm1, 58, 59)
    set_conv_weight(model.layer5.conv1, 60, 3, 96, 96)
    set_norm_weight(model.layer5.norm2, 63, 64)
    set_conv_weight(model.layer5.conv2, 65, 3, 96, 96)

    #final_norm
    set_norm_weight(model.final_norm, 67, 68)
    set_conv_weight(model.p_conv1, 69, 1, 96, 32)
    set_conv_weight(model.p_conv2, 70, 1, 96, 32)
    set_norm_weight(model.p_norm2, 72, 73)
    set_linear_weight(model.p_linear_g, 74, 96, 32)
    set_norm_weight(model.p_norm_combine, 76, 77)
    set_conv_weight(model.p_conv_final, 78, 1, 32, 2)


model = KataNet()
load_weights_from_bins(model, bins)
model.eval() #必须要加这一行否则diff直接几千
output = model(input19, inputGlobal19, 19)
diff(output, policy19)


input9=torch.tensor(input9, dtype=torch.float32).unsqueeze(0)
inputGlobal9=torch.tensor(inputGlobal9, dtype=torch.float32).unsqueeze(0)
output = model(input9, inputGlobal9, 9)

policy9=torch.tensor(policy9, dtype=torch.float32).unsqueeze(0)
diff(output, policy9)