#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>

typedef struct {
    char* content;
    int text_start;
    int text_len;
    int bin_start;
    int bin_len;
    float *floats;
} bin_t;


bin_t BINS[]={
NULL, 0, 169, 169, 0, NULL,
NULL, 76201, 32, 76233, 0, NULL,
NULL, 83529, 82, 83611, 0, NULL,
NULL, 83995, 6, 84001, 0, NULL,
NULL, 84385, 6, 84391, 0, NULL,
NULL, 84775, 101, 84876, 0, NULL,
NULL, 416652, 52, 416704, 0, NULL,
NULL, 417088, 6, 417094, 0, NULL,
NULL, 417478, 6, 417484, 0, NULL,
NULL, 417868, 6, 417874, 0, NULL,
NULL, 418258, 101, 418359, 0, NULL,
NULL, 750135, 82, 750217, 0, NULL,
NULL, 750601, 6, 750607, 0, NULL,
NULL, 750991, 6, 750997, 0, NULL,
NULL, 751381, 101, 751482, 0, NULL,
NULL, 1083258, 52, 1083310, 0, NULL,
NULL, 1083694, 6, 1083700, 0, NULL,
NULL, 1084084, 6, 1084090, 0, NULL,
NULL, 1084474, 6, 1084480, 0, NULL,
NULL, 1084864, 101, 1084965, 0, NULL,
NULL, 1416741, 79, 1416820, 0, NULL,
NULL, 1417204, 6, 1417210, 0, NULL,
NULL, 1417594, 6, 1417600, 0, NULL,
NULL, 1417984, 112, 1418096, 0, NULL,
NULL, 1639280, 64, 1639344, 0, NULL,
NULL, 1749936, 62, 1749998, 0, NULL,
NULL, 1750126, 6, 1750132, 0, NULL,
NULL, 1750260, 6, 1750266, 0, NULL,
NULL, 1750394, 116, 1750510, 0, NULL,
NULL, 1775086, 52, 1775138, 0, NULL,
NULL, 1775394, 6, 1775400, 0, NULL,
NULL, 1775656, 6, 1775662, 0, NULL,
NULL, 1775918, 6, 1775924, 0, NULL,
NULL, 1776180, 101, 1776281, 0, NULL,
NULL, 1997465, 82, 1997547, 0, NULL,
NULL, 1997931, 6, 1997937, 0, NULL,
NULL, 1998321, 6, 1998327, 0, NULL,
NULL, 1998711, 101, 1998812, 0, NULL,
NULL, 2330588, 52, 2330640, 0, NULL,
NULL, 2331024, 6, 2331030, 0, NULL,
NULL, 2331414, 6, 2331420, 0, NULL,
NULL, 2331804, 6, 2331810, 0, NULL,
NULL, 2332194, 101, 2332295, 0, NULL,
NULL, 2664071, 79, 2664150, 0, NULL,
NULL, 2664534, 6, 2664540, 0, NULL,
NULL, 2664924, 6, 2664930, 0, NULL,
NULL, 2665314, 112, 2665426, 0, NULL,
NULL, 2886610, 64, 2886674, 0, NULL,
NULL, 2997266, 62, 2997328, 0, NULL,
NULL, 2997456, 6, 2997462, 0, NULL,
NULL, 2997590, 6, 2997596, 0, NULL,
NULL, 2997724, 116, 2997840, 0, NULL,
NULL, 3022416, 52, 3022468, 0, NULL,
NULL, 3022724, 6, 3022730, 0, NULL,
NULL, 3022986, 6, 3022992, 0, NULL,
NULL, 3023248, 6, 3023254, 0, NULL,
NULL, 3023510, 101, 3023611, 0, NULL,
NULL, 3244795, 82, 3244877, 0, NULL,
NULL, 3245261, 6, 3245267, 0, NULL,
NULL, 3245651, 6, 3245657, 0, NULL,
NULL, 3246041, 101, 3246142, 0, NULL,
NULL, 3577918, 52, 3577970, 0, NULL,
NULL, 3578354, 6, 3578360, 0, NULL,
NULL, 3578744, 6, 3578750, 0, NULL,
NULL, 3579134, 6, 3579140, 0, NULL,
NULL, 3579524, 101, 3579625, 0, NULL,
NULL, 3911401, 41, 3911442, 0, NULL,
NULL, 3911826, 6, 3911832, 0, NULL,
NULL, 3912216, 6, 3912222, 0, NULL,
NULL, 3912606, 100, 3912706, 0, NULL,
NULL, 3924994, 45, 3925039, 0, NULL,
NULL, 3937327, 43, 3937370, 0, NULL,
NULL, 3937498, 6, 3937504, 0, NULL,
NULL, 3937632, 6, 3937638, 0, NULL,
NULL, 3937766, 78, 3937844, 0, NULL,
NULL, 3950132, 43, 3950175, 0, NULL,
NULL, 3950303, 6, 3950309, 0, NULL,
NULL, 3950437, 6, 3950443, 0, NULL,
NULL, 3950571, 83, 3950654, 0, NULL,
NULL, 3950910, 42, 3950952, 0, NULL,
NULL, 3963240, 44, 3963284, 0, NULL,
NULL, 3963412, 85, 3963497, 0, NULL,
NULL, 3963753, 60, 3963813, 0, NULL,
NULL, 3976101, 42, 3976143, 0, NULL,
NULL, 3976271, 6, 3976277, 0, NULL,
NULL, 3976405, 6, 3976411, 0, NULL,
NULL, 3976539, 75, 3976614, 0, NULL,
NULL, 4001190, 32, 4001222, 0, NULL,
NULL, 4001478, 83, 4001561, 0, NULL,
NULL, 4002329, 40, 4002369, 0, NULL,
NULL, 4002381, 49, 4002430, 0, NULL,
NULL, 4003966, 44, 4004010, 0, NULL,
NULL, 4004034, 51, 4004085, 0, NULL

};

extern const float trunkScratch[96][19][19];
extern const float input[22][19][19];
extern const float inputGlobal[19];
extern const float inputMatMulOut[96];
extern const float trunkScratch_afterBias[96][19][19];
extern const float afternorm[96][19][19];
extern const float afterconv[96][19][19];
extern const float afternorm2[96][19][19];
extern const float afterconv2[96][19][19];
extern const float after2ordis[96][19][19];
extern const float regularOut[64][19][19];
extern const float gpoolOut2[32][19][19];
extern const float afterglobal[96][19][19];
extern const float afterblocks[96][19][19];
extern const float policy[2][19][19];

long get_file_size(FILE *fp) {
    long cur = ftell(fp);          // 记录当前位置
    fseek(fp, 0, SEEK_END);        // 跳到文件末尾
    long size = ftell(fp);         // 获得文件大小（字节）
    fseek(fp, cur, SEEK_SET);      // 恢复文件指针位置
    return size;
}

char *read_ascii_name(FILE *fp, int offset, int len) {
    unsigned char *buf = malloc(len + 1);
    if (!buf) return NULL;

    fseek(fp, offset, SEEK_SET);
    fread(buf, 1, len, fp);
    buf[len] = 0;

    // 删除不可打印字符，例如 '\n', 0x00, < > ? 等
    for (int i = 0; i < len; i++) {
        if (buf[i] < 32 || buf[i] > 126)
            buf[i] = '|';  // 用 0 截断
    }

    buf[len-5] = 0; // @BIN@不输出

    return (char*)buf;
}

int float_cmp(const void* a, const void* b) {
    float fa = *(const float*)a;
    float fb = *(const float*)b;
    if (fa < fb) return -1;
    else if (fa > fb) return 1;
    else return 0;
}

float focus22(float input[22][3][3], float kernel[22][3][3])
{
  float s = 0.0f;
  for(int i =0; i<22; i++)
  {
    for(int j=0; j<3; j++)
    {
      for(int k=0; k<3; k++)
      {
        if(input[i][j][k] > 1000000.0f || kernel[i][j][k] > 1000000.0f)
        {
          printf("warning! big number");
        }
        s += input[i][j][k] * kernel[i][j][k];
      }
    }
  }
  return s;
}

void fries22(float input_padded[22][21][21], float output[22][3][3], int x, int y)
{
  for(int i=0; i<22; i++)
  {
    for(int j=0; j<3; j++)
    {
      for(int k=0; k<3; k++)
      {
        output[i][j][k] = input_padded[i][x+j][y+k];
        if(output[i][j][k] > 1000000.0f)
        {
          printf("warning! big number");
        }
      }
    }
  }

}

void slip22(const float input[22][19][19], float output[19][19], float kernel[22][3][3])
{
  // padding
  float input_padded[22][19+1+1][19+1+1] = {0};
  for(int i=0; i<22; i++)
  {
    for(int j=0; j<19; j++)
    {
      for(int k=0; k<19; k++)
      {
        input_padded[i][j+1][k+1] = input[i][j][k];
      }
    }
  }
  // conv
  for(int i=0; i<19; i++)
  {
    for(int j=0; j<19; j++)
    {
      float fry[22][3][3];
      fries22(input_padded, fry, i, j);
      output[i][j] = focus22(fry, kernel);
    }
  }
}

void conv2296(const float input[22][19][19], float output[96][19][19], float kernel[96][22][3][3])
{
  for (int I=0; I<96; I++)
  {
    float (*cannon)[3][3] = kernel[I]; /*cannon[22][3][3]*/
    slip22(input, output[I], cannon);
  }
}

float focus96(float input[96][3][3], float kernel[96][3][3])
{
  float s = 0.0f;
  for(int i =0; i<96; i++)
  {
    for(int j=0; j<3; j++)
    {
      for(int k=0; k<3; k++)
      {
        if(input[i][j][k] > 1000000.0f || kernel[i][j][k] > 1000000.0f)
        {
          printf("warning! big number");
        }
        s += input[i][j][k] * kernel[i][j][k];
      }
    }
  }
  return s;
}

void fries96(float input_padded[96][21][21], float output[96][3][3], int x, int y)
{
  for(int i=0; i<96; i++)
  {
    for(int j=0; j<3; j++)
    {
      for(int k=0; k<3; k++)
      {
        output[i][j][k] = input_padded[i][x+j][y+k];
        if(output[i][j][k] > 1000000.0f)
        {
          printf("warning! big number");
        }
      }
    }
  }

}

void slip96(const float input[96][19][19], float output[19][19], float kernel[96][3][3])
{
  // padding
  float input_padded[96][19+1+1][19+1+1] = {0};
  for(int i=0; i<96; i++)
  {
    for(int j=0; j<19; j++)
    {
      for(int k=0; k<19; k++)
      {
        input_padded[i][j+1][k+1] = input[i][j][k];
      }
    }
  }
  // conv
  for(int i=0; i<19; i++)
  {
    for(int j=0; j<19; j++)
    {
      float fry[96][3][3];
      fries96(input_padded, fry, i, j);
      output[i][j] = focus96(fry, kernel);
    }
  }
}

float focus64(float input[64][3][3], float kernel[64][3][3])
{
  float s = 0.0f;
  for(int i =0; i<64; i++)
  {
    for(int j=0; j<3; j++)
    {
      for(int k=0; k<3; k++)
      {
        if(input[i][j][k] > 1000000.0f || kernel[i][j][k] > 1000000.0f)
        {
          printf("warning! big number");
        }
        s += input[i][j][k] * kernel[i][j][k];
      }
    }
  }
  return s;
}

void fries64(float input_padded[64][21][21], float output[64][3][3], int x, int y)
{
  for(int i=0; i<64; i++)
  {
    for(int j=0; j<3; j++)
    {
      for(int k=0; k<3; k++)
      {
        output[i][j][k] = input_padded[i][x+j][y+k];
        if(output[i][j][k] > 1000000.0f)
        {
          printf("warning! big number");
        }
      }
    }
  }

}

void slip64(const float input[64][19][19], float output[19][19], float kernel[64][3][3])
{
  // padding
  float input_padded[64][19+1+1][19+1+1] = {0};
  for(int i=0; i<64; i++)
  {
    for(int j=0; j<19; j++)
    {
      for(int k=0; k<19; k++)
      {
        input_padded[i][j+1][k+1] = input[i][j][k];
      }
    }
  }
  // conv
  for(int i=0; i<19; i++)
  {
    for(int j=0; j<19; j++)
    {
      float fry[64][3][3];
      fries64(input_padded, fry, i, j);
      output[i][j] = focus64(fry, kernel);
    }
  }
}

void conv9696(const float input[96][19][19], float output[96][19][19], float kernel[96][96][3][3])
{
  for (int I=0; I<96; I++)
  {
    float (*cannon)[3][3] = kernel[I]; /*cannon[22][3][3]*/
    slip96(input, output[I], cannon);
  }
}

void conv9664(const float input[96][19][19], float output[64][19][19], float kernel[64][96][3][3])
{
  for (int I=0; I<64; I++)
  {
    float (*cannon)[3][3] = kernel[I]; /*cannon[22][3][3]*/
    slip96(input, output[I], cannon);
  }
}

void conv9632(const float input[96][19][19], float output[32][19][19], float kernel[32][96][3][3])
{
  for (int I=0; I<32; I++)
  {
    float (*cannon)[3][3] = kernel[I]; /*cannon[22][3][3]*/
    slip96(input, output[I], cannon);
  }
}

void conv6496(const float input[64][19][19], float output[96][19][19], float kernel[96][64][3][3])
{
  for (int I=0; I<96; I++)
  {
    float (*cannon)[3][3] = kernel[I]; /*cannon[22][3][3]*/
    slip64(input, output[I], cannon);
  }
}

void conv1x1(const float input[96][19][19], float output[32][19][19], float kernel[32][96])
{
  for (int I=0; I<32; I++)
  {
    for(int i=0; i<19; i++)
      for(int j=0; j<19; j++)
      {
        float s = 0.0f;
        for (int k=0; k<96; k++)
        {
          s+=input[k][i][j]*kernel[I][k];
        }
        output[I][i][j] = s;
      }
  }
}

void conv322(const float input[32][19][19], float output[2][19][19], float kernel[2][32])
{
  for (int I=0; I<2; I++)
  {
    for(int i=0; i<19; i++)
      for(int j=0; j<19; j++)
      {
        float s = 0.0f;
        for (int k=0; k<32; k++)
        {
          s+=input[k][i][j]*kernel[I][k];
        }
        output[I][i][j] = s;
      }
  }
}

float err1(float *mat1, const float *mat2, int x, float *maxdiff, int *oi)
{
  float s = 0.0f;
  float max_diff = -999.0f;
  for (int i=0; i<x; i++)
  {
    float a = mat1[i] - mat2[i];
    float diff;
    if(a > 0)
      diff = a;
    else
      diff = -a;
  
    if(diff > 0.1)
    {
      printf("error: diff too big!! %d: %f\n", i, diff);
      exit(EXIT_FAILURE);
    }
  
    s += diff;
    if(diff > max_diff)
    {
      max_diff=diff;
      *oi = i;
    }
    if(s>1000000.0f)
    {
      printf("warning! big number");
    }
  }
  *maxdiff = max_diff;
  return s;
}

float err3(float *mat1, const float *mat2, int x, int y, int z, float *maxdiff, int *oi, int *oj, int *ok)
{
  float s = 0.0f;
  float max_diff = -999.0f;
  for (int i=0; i<x; i++)
    for (int j=0; j<y; j++)
      for (int k=0; k<z; k++)
      {
        int idx = i * y * z + j * z + k;
        float a = mat1[idx] - mat2[idx];
        float diff;
        if(a > 0)
          diff = a;
        else
          diff = -a;

        if(diff > 0.001)
        {
          printf("error: diff too big!! %d, (%d, %d): output: %f, standard: %f, diff: %f\n", i, j, k, mat1[idx], mat2[idx], diff);
          exit(EXIT_FAILURE);
        }

        s += diff;
        if(diff > max_diff)
        {
          max_diff=diff;
          *oi = i;
          *oj = j;
          *ok = k;
        }
        if(s>1000000.0f)
        {
          printf("warning! big number");
        }
      }
  *maxdiff = max_diff;
  return s;
}

void calc1(const float input[19], float output[96], float linear[96][19])
{
  for (int i = 0; i < 96; i++) {
    float sum = 0.0f;
    for (int j = 0; j < 19; j++) {
        sum += linear[i][j] * input[j];
    }
    output[i] = sum;
  }
}

void calcBias(const float input[96][19][19], float output[96][19][19], float bias[96])
{
  for (int i = 0; i < 96; i++) {
    for (int j = 0; j < 19; j++) 
      for (int k = 0; k < 19; k++){

          output[i][j][k] = input[i][j][k] + bias[i];
      }
  }
}

void norm(const float input[96][19][19], float output[96][19][19], float scale[96], float bias[96])
{
    for (int i = 0; i < 96; i++) {
      for (int j = 0; j < 19; ++j) {
        for (int k = 0; k < 19; ++k) {
          float x = input[i][j][k] * scale[i] + bias[i];
          output[i][j][k] = (x > 0.0f) ? x : 0.0f;

          // 小棋盘注意过滤掉不在棋盘上的点
          }
        }
    }
}

void norm32(const float input[32][19][19], float output[32][19][19], float scale[32], float bias[32])
{
    for (int i = 0; i < 32; i++) {
      for (int j = 0; j < 19; ++j) {
        for (int k = 0; k < 19; ++k) {
          float x = input[i][j][k] * scale[i] + bias[i];
          output[i][j][k] = (x > 0.0f) ? x : 0.0f;

          // 小棋盘注意过滤掉不在棋盘上的点
          }
        }
    }
}

void norm64(const float input[64][19][19], float output[64][19][19], float scale[64], float bias[64])
{
    for (int i = 0; i < 64; i++) {
      for (int j = 0; j < 19; ++j) {
        for (int k = 0; k < 19; ++k) {
          float x = input[i][j][k] * scale[i] + bias[i];
          output[i][j][k] = (x > 0.0f) ? x : 0.0f;

          // 小棋盘注意过滤掉不在棋盘上的点
          }
        }
    }
}

void add(const float input[96][19][19], float output[96][19][19], float adder[96][19][19])
{
  for (int i = 0; i < 96; i++) {
    for (int j = 0; j < 19; ++j) {
      for (int k = 0; k < 19; ++k) {
          output[i][j][k] = input[i][j][k] + adder[i][j][k];
        }
      }
  }
}

void add_broadcast(const float input[64][19][19], float output[64][19][19], float adder[64])
{
  for (int i = 0; i < 64; i++) {
    for (int j = 0; j < 19; ++j) {
      for (int k = 0; k < 19; ++k) {
        output[i][j][k] = input[i][j][k] + adder[i];
        }
      }
  }
}

void add_broadcast32(const float input[32][19][19], float output[32][19][19], float adder[32])
{
  for (int i = 0; i < 32; i++) {
    for (int j = 0; j < 19; ++j) {
      for (int k = 0; k < 19; ++k) {
        output[i][j][k] = input[i][j][k] + adder[i];
        }
      }
  }
}

void rowsG(float input[32][19][19], float output[96])
{
  float div = 361.0f;
  float sqrtdiv = 19.0f;
  for(int I = 0; I<32; I++)
  {
    float s = 0.0f;
    float max = -1.0f;
    for(int i=0; i<19; i++)
      for(int j=0; j<19; j++)
      {
        float x = input[I][i][j];
        s += x;
        /* remember change this*/
        float maskVal = 1.0f;
        float temp = x + (maskVal - 1.0f);
        max = temp > max ? temp : max;
      }
        
    float mean = s / div;
    output[0 + I] = mean;
    output[32 + I] = mean * (sqrtdiv - 14.0f) * 0.1f;
    output[64 + I] = max;
  }

}

void linear(const float input[96], float output[64], float nn[64][96])
{
    for (int I = 0; I < 64; I++) {
      float s = 0.0f;
      for (int i = 0; i < 96; i++) {
          s += input[i] * nn[I][i];
      }
      output[I] = s;
    }
}

void linear9632(const float input[96], float output[32], float nn[32][96])
{
    for (int I = 0; I < 32; I++) {
      float s = 0.0f;
      for (int i = 0; i < 96; i++) {
          s += input[i] * nn[I][i];
      }
      output[I] = s;
    }
}

void ordi(float input[96][19][19], float output[96][19][19], float scale0[96], float bias0[96], const float (*afternorm)[19][19], float kernel1[96][96][3][3], const float (*afterconv)[19][19], float scale1[96], float bias1[96], const float (*afternorm2)[19][19], float kernel2[96][96][3][3], const float (*afterconv2)[19][19])
{
    float maxdiff;
    int erri, errj, errk;

    float output0[96][19][19];

    norm(input, output0, scale0, bias0);

    if(afternorm != NULL)
    {
      printf("sumdiff: %f\n", err3((float*)output0, (const float*)afternorm, 96, 19, 19, &maxdiff, &erri, &errj, &errk));
      printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);
    }
    
    float output1[96][19][19];

    conv9696(output0, output1, kernel1);

    if(afterconv != NULL)
    {
      printf("sumdiff: %f\n", err3((float*)output1, (const float*)afterconv, 96, 19, 19, &maxdiff, &erri, &errj, &errk));
      printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);
    }
    
    float output2[96][19][19];

    norm(output1, output2, scale1, bias1);

    if(afternorm2 != NULL)
    {
      printf("sumdiff: %f\n", err3((float*)output2, (const float*)afternorm2, 96, 19, 19, &maxdiff, &erri, &errj, &errk));
      printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);
    }
    float output3[96][19][19];

    conv9696(output2, output3, kernel2);

    add(input, output, output3);

    if(afterconv2 != NULL)
    {
      printf("sumdiff: %f\n", err3((float*)output, (const float*)afterconv2, 96, 19, 19, &maxdiff, &erri, &errj, &errk));
      printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);
    }

}

void gpool(float input[96][19][19], float output[96][19][19], float scale0[96], float bias0[96], float kernel0[64][96][3][3], const float (*regularOut)[19][19], float kernel1[32][96][3][3], float scale1[32], float bias1[32], const float (*gpoolOut)[19][19], float nn[64][96], float scale2[64], float bias2[64], float kernel2[96][64][3][3], const float (*afterglobal)[19][19])
{
    float maxdiff;
    int erri, errj, errk;

    float output0[96][19][19];

    norm(input, output0, scale0, bias0);

    //从这开始分为了两支，先是regular支

    float output1[64][19][19];

    conv9664(output0, output1, kernel0);

    if(regularOut != NULL)
    {
      printf("sumdiff: %f\n", err3((float*)output1, (const float*)regularOut, 64, 19, 19, &maxdiff, &erri, &errj, &errk));
      printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);
    }

    // g分支

    float output2[32][19][19];

    conv9632(output0, output2, kernel1);

    float output3[32][19][19];

    norm32(output2, output3, scale1, bias1);

    if(gpoolOut != NULL)
    {
      printf("sumdiff: %f\n", err3((float*)output3, (const float*)gpoolOut, 32, 19, 19, &maxdiff, &erri, &errj, &errk));
      printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);
    }
    
    float output4[96];
    rowsG(output3, output4);

    /* 1.030514, 0.515257, 7.487195 */
    printf("output4: %f, %f, %f\n", output4[0], output4[32], output4[64]);

    float output5[64];
    linear(output4, output5, nn);

    /* -2.681623, 0.782912 */
    printf("output5: %f, %f\n", output5[0], output5[55]);


    float output6[64][19][19];
    add_broadcast(output1, output6, output5);

    
    // 汇聚后再来一次normconv
    // BINS[29]全是0，跳过
    // BINS[30]跳过

    float output7[64][19][19];

    norm64(output6, output7, scale2, bias2);

    float output8[96][19][19];

    conv6496(output7, output8, kernel2);

    add(input, output, output8);

    if(afterglobal != NULL)
    {
      printf("sumdiff: %f\n", err3((float*)output, (const float*)afterglobal, 96, 19, 19, &maxdiff, &erri, &errj, &errk));
      printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);
    }
}

void load_ordi(float scale0[96], float bias0[96], float kernel1[96][96][3][3], float scale1[96], float bias1[96], float kernel2[96][96][3][3], int i0, int i1, int i2, int i3, int i4, int i5)
{
  int n=0;

  for(int i=0; i < 96; i++)
  {
    scale0[i] = BINS[i0].floats[n];
    bias0[i] = BINS[i1].floats[n];
    n++;
  }

  n=0;

  for(int k=0; k <3; k++)
    for(int l=0; l <3; l++)
      for(int j=0; j <96; j++)
        for(int i=0; i <96; i++)
        {
          kernel1[i][j][k][l] = BINS[i2].floats[n++];
        }
  
  n=0;

  for(int i=0; i < 96; i++)
  {
    scale1[i] = BINS[i3].floats[n];
    bias1[i] = BINS[i4].floats[n];
    n++;
  }

  n=0;

  for(int k=0; k <3; k++)
    for(int l=0; l <3; l++)
      for(int j=0; j <96; j++)
        for(int i=0; i <96; i++)
        {
          kernel2[i][j][k][l] = BINS[i5].floats[n++];
        }
}

void load_gpool(float scale0[96], float bias0[96], float kernel0[64][96][3][3], float kernel1[32][96][3][3], float scale1[32], float bias1[32], float nn[64][96], float scale2[64], float bias2[64], float kernel2[96][64][3][3], int i0, int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9)
{
    int n=0;

    for(int i=0; i < 96; i++)
    {
      scale0[i] = BINS[i0].floats[n];
      bias0[i] = BINS[i1].floats[n];
      n++;
    }

    n=0;

    for(int k=0; k <3; k++)
      for(int l=0; l <3; l++)
        for(int j=0; j <96; j++)
          for(int i=0; i <64; i++)
          {
            kernel0[i][j][k][l] = BINS[i2].floats[n++];
          }

    n=0;

    for(int k=0; k <3; k++)
      for(int l=0; l <3; l++)
        for(int j=0; j <96; j++)
          for(int i=0; i <32; i++)
          {
            kernel1[i][j][k][l] = BINS[i3].floats[n++];
          }

    n=0;

    for(int i=0; i < 32; i++)
    {
      scale1[i] = BINS[i4].floats[n];
      bias1[i] = BINS[i5].floats[n];
      n++;
    }

    n=0;

    for(int j=0; j < 96; j++)
      for(int i=0; i < 64; i++)
      {
        nn[i][j] = BINS[i6].floats[n];
        n++;
      }

    n=0;

    for(int i=0; i < 64; i++)
    {
      scale2[i] = BINS[i7].floats[n];
      bias2[i] = BINS[i8].floats[n];
      n++;
    }

    n=0;

    for(int k=0; k <3; k++)
      for(int l=0; l <3; l++)
        for(int j=0; j <64; j++)
          for(int i=0; i <96; i++)
          {
            kernel2[i][j][k][l] = BINS[i9].floats[n++];
          }
}


int main() {
    const char *gzfile = "model3e4.bin.gz";
    const char *binfile = "model3e4.bin";

    /* 1. 打开 gzip 文件 */
    gzFile in = gzopen(gzfile, "rb");
    if (!in) {
        fprintf(stderr, "无法打开 %s\n", gzfile);
        return 1;
    }

    /* 2. 打开输出文件 */
    FILE *out = fopen(binfile, "wb");
    if (!out) {
        fprintf(stderr, "无法创建 %s\n", binfile);
        gzclose(in);
        return 1;
    }

    /* 3. 解压写入 model3e4.bin */
    char buffer[4096];
    int bytes;
    while ((bytes = gzread(in, buffer, sizeof(buffer))) > 0) {
        fwrite(buffer, 1, bytes, out);
    }

    gzclose(in);
    fclose(out);

    printf("unzip finished: %s\n", binfile);

    // 4. 打开解压出的 .bin 文件，按 float32 读取 
    FILE *fp = fopen(binfile, "rb");
    if (!fp) {
        fprintf(stderr, "无法打开 %s\n", binfile);
        return 1;
    }

    long size = get_file_size(fp);

    int N = sizeof(BINS)/sizeof(bin_t);

    for (int i = 0; i < N; i++) {
        bin_t bin = BINS[i];
        if(bin.text_start + bin.text_len != bin.bin_start)
        {
            fprintf(stderr, "text len check wrong, doesn't match bin start: %d + %d != %d\n", bin.text_start, bin.text_len, bin.bin_start);
            return 1;
        }
    }

    printf("text len check success, matched bin start.\n");

    for (int i = 0; i < N; i++) {
        if(i != N-1)
        {

            BINS[i].bin_len = BINS[i+1].text_start - BINS[i].bin_start;
        }
        else
        {
            BINS[i].bin_len = (size-1) /*末尾有一个0xa*/ - BINS[i].bin_start;
        }

        if(BINS[i].bin_len %4 != 0)
        {
            fprintf(stderr, "bin len check wrong, doesn't 4x: %d, %d, %d, %d\n", BINS[i].text_start, BINS[i].text_len, BINS[i].bin_start, BINS[i].bin_len);
            return 1;
        }
    }

    printf("bin len check success, matched 4x.\n");

    //填充content

    for (int i = 0; i < N; i++) {
        BINS[i].content = read_ascii_name(fp, BINS[i].text_start, BINS[i].text_len);
        if (!BINS[i].content) {
            fprintf(stderr, "fail to read name for %d\n", i);
        }
    }

    size_t total_floats = 0;
    float total_min = 999;
    float total_max = -999;

    /* 按 BINS 区间逐段读取 float */
    for (int i = 0; i < N; i++) {
        bin_t *b = &BINS[i];

        int float_count = b->bin_len / 4;
        float *data = malloc(sizeof(float) * float_count);
        if(!data) { fprintf(stderr, "内存分配失败\n"); return 1; }

        /* 读取该区间 */
        fseek(fp, b->bin_start, SEEK_SET);
        fread(data, sizeof(float), float_count, fp);

        // 计算平均值
        double sum = 0;
        float min = data[0];
        float max = data[0];
        b->floats = malloc(sizeof(float) * float_count);
        for(int j=0; j<float_count; j++) {
            float v = data[j];
            sum += v;
            if(v < min) min = v;
            if(v > max) max = v;
            b->floats[j] = v;
        }
        double mean = sum / float_count;

        // 计算中位数
        qsort(data, float_count, sizeof(float), float_cmp);
        double median = (float_count % 2 == 0) ?
            (data[float_count/2 - 1] + data[float_count/2])/2.0 :
            data[float_count/2];

        printf("BIN %d: content=%s, floats=%d, min=%.6f, max=%.6f, mean=%.6f, median=%.6f\n",
                i, b->content, float_count, min, max, mean, median);
        free(data);

        total_floats += float_count;
        if(min < total_min) total_min = min;
        if(max > total_max) total_max = max;
    }

    printf("%d\n", BINS[0].bin_len / 4);

    printf("========\ntotal float amount: %zu, total min:%.6f, total max:%.6f\n", total_floats, total_min, total_max);
    fclose(fp);

    float output0[96][19][19];
    float kernel0[96][22][3][3];
    int n = 0;
 
    for(int k=0; k <3; k++)
      for(int l=0; l <3; l++)
        for(int j=0; j <22; j++)
          for(int i=0; i <96; i++)
          {
            kernel0[i][j][k][l] = BINS[0].floats[n++];
          }
    conv2296(input, output0, kernel0);

    float maxdiff;
    int erri, errj, errk;

    printf("sumdiff: %f\n", err3((float*)output0, (const float*)trunkScratch, 96, 19, 19, &maxdiff, &erri, &errj, &errk));

    printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);



    float output1[96];
    float linear0[96][19];
    n=0;

    for(int j=0; j <19; j++)
      for(int i=0; i < 96; i++)
      {
        linear0[i][j] = BINS[1].floats[n++];
      }

    calc1(inputGlobal, output1, linear0);

    
    printf("sumdiff: %f\n", err1((float*)output1, (const float*)inputMatMulOut, 96, &maxdiff, &erri));

    printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);


    float output2[96][19][19];
    calcBias(output0, output2, output1);

    printf("sumdiff: %f\n", err3((float*)output2, (const float*)trunkScratch_afterBias, 96, 19, 19, &maxdiff, &erri, &errj, &errk));

    printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);

    // block阶段
    // 6个block分别是 ordi, ordi, gpool, ordi, gpool, ordi

    float scale0[96],            scale2[96],            scale7[96],            scale12[96];
    float bias0[96],             bias2[96],             bias7[96],             bias12[96];
    float kernel1[96][96][3][3], kernel3[96][96][3][3], kernel8[96][96][3][3],  kernel13[96][96][3][3];
    float scale1[96],            scale3[96],            scale8[96],            scale13[96];
    float bias1[96],             bias3[96],             bias8[96],             bias13[96];
    float kernel2[96][96][3][3], kernel4[96][96][3][3], kernel9[96][96][3][3], kernel14[96][96][3][3];

    float scale4[96],            scale9[96];
    float bias4[96],             bias9[96];
    float kernel5[64][96][3][3], kernel10[64][96][3][3];
    float kernel6[32][96][3][3], kernel11[32][96][3][3];
    float scale5[32],            scale10[32];
    float bias5[32],             bias10[32];
    float nn0[64][96],           nn1[64][96];
    float scale6[64],            scale11[64];
    float bias6[64],             bias11[64];
    float kernel7[96][64][3][3], kernel12[96][64][3][3];

    float output3[96][19][19];
    float output4[96][19][19];
    float output5[96][19][19];
    float output6[96][19][19];
    float output7[96][19][19];
    float output8[96][19][19];

    load_ordi(scale0, bias0, kernel1, scale1, bias1, kernel2, 3, 4, 5, 8, 9, 10);
    load_ordi(scale2, bias2, kernel3, scale3, bias3, kernel4, 12, 13, 14, 17, 18, 19);
    load_gpool(scale4, bias4, kernel5, kernel6, scale5, bias5, nn0, scale6, bias6, kernel7, 21, 22, 23, 24, 26, 27, 28, 31, 32, 33);
    load_ordi(scale7, bias7, kernel8, scale8, bias8, kernel9, 35, 36, 37, 40, 41, 42);
    load_gpool(scale9, bias9, kernel10, kernel11, scale10, bias10, nn1, scale11, bias11, kernel12, 44, 45, 46, 47, 49, 50, 51, 54, 55, 56);
    load_ordi(scale12, bias12, kernel13, scale13, bias13, kernel14, 58, 59, 60, 63, 64, 65);

    ordi (output2, output3, scale0,  bias0,  afternorm, kernel1,    afterconv, scale1,  bias1,  afternorm2, kernel2,  afterconv2);
    ordi (output3, output4, scale2,  bias2,  NULL,      kernel3,    NULL,      scale3,  bias3,  NULL,       kernel4,  after2ordis);
    gpool(output4, output5, scale4,  bias4,  kernel5,   regularOut, kernel6,   scale5,  bias5,  gpoolOut2,  nn0,      scale6,  bias6,  kernel7,  afterglobal);
    ordi (output5, output6, scale7,  bias7,  NULL,      kernel8,    NULL,      scale8,  bias8,  NULL,       kernel9,  NULL);
    gpool(output6, output7, scale9,  bias9,  kernel10,  NULL,       kernel11,  scale10, bias10, NULL,       nn1,      scale11, bias11, kernel12, NULL);
    ordi (output7, output8, scale12, bias12, NULL,      kernel13,   NULL,      scale13, bias13, NULL,       kernel14, NULL);

    printf("sumdiff: %f\n", err3((float*)output8, (const float*)afterblocks, 96, 19, 19, &maxdiff, &erri, &errj, &errk));
    printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);

    /* 6个block全部结束*/


    float scale14[96];
    float bias14[96];
    float output9[96][19][19];
    n=0;

    for(int i=0; i < 96; i++)
    {
      scale14[i] = BINS[67].floats[n];
      bias14[i] = BINS[68].floats[n];
      n++;
    }

    norm(output8, output9, scale14, bias14);

    /* policy头 */

    /*分支1备用*/
    float kernel15[32][96];

    n=0;

    for(int j=0; j <96; j++)
      for(int i=0; i <32; i++)
      {
        kernel15[i][j] = BINS[69].floats[n++];
      }

    float output10[32][19][19];
    conv1x1(output9, output10, kernel15);

    /*分支2*/

    float kernel16[32][96];

    n=0;

    for(int j=0; j <96; j++)
      for(int i=0; i <32; i++)
      {
        kernel16[i][j] = BINS[70].floats[n++];
      }

    float output11[32][19][19];
    conv1x1(output9, output11, kernel16);

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

    norm32(output11, output12, scale15, bias15);

    float output13[96];
    rowsG(output12, output13);

    float nn2[32][96];
    n=0;

    for(int j=0; j < 96; j++)
      for(int i=0; i < 32; i++)
      {
        nn2[i][j] = BINS[74].floats[n];
        n++;
      }

    float output14[32];
    linear9632(output13, output14, nn2);

    float output15[32][19][19];
    add_broadcast32(output10, output15, output14);

    /* output15[0][0][0]: 0.358652, output15[0][0][1]: -0.280337, output15[0][0][2]: 0.398067 */
    printf("output15[0][0][0]: %f, output15[0][0][1]: %f, output15[0][0][2]: %f\n", output15[0][0][0], output15[0][0][1], output15[0][0][2]);

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

    norm32(output15, output16, scale16, bias16);

    float kernel17[2][32];

    n=0;

    for(int j=0; j <32; j++)
      for(int i=0; i <2; i++)
      {
        kernel17[i][j] = BINS[78].floats[n++];
      }
    
    float output17[2][19][19];

    conv322(output16, output17, kernel17);

    /* 4.33436871 */
    printf("output17[0][1][14]: %f\n", output17[0][1][14]);

    printf("sumdiff: %f\n", err3((float*)output17, (const float*)policy, 2, 19, 19, &maxdiff, &erri, &errj, &errk));
    printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);


    














    return 0;
}
