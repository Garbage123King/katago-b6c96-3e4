#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>

#define TRUE 1
#define FALSE 0

typedef struct {
    char* content;
    int text_start;
    int text_len;
    int bin_start;
    int bin_len;
    float *floats;
} bin_t;

typedef struct
{
  float scale0[96];
  float bias0[96];
  float kernel0[96][96][3][3];
  float scale1[96];
  float bias1[96];
  float kernel1[96][96][3][3];
} ordi_param;

typedef struct
{
  float scale0[96];
  float bias0[96];
  float kernel0[64][96][3][3];
  float kernel1[32][96][3][3];
  float scale1[32];
  float bias1[32];
  float nn0[64][96];
  float scale2[64];
  float bias2[64];
  float kernel2[96][64][3][3];
} gpool_param;

typedef struct
{
  ordi_param block0;
  ordi_param block1;
  gpool_param block2;
  ordi_param block3;
  gpool_param block4;
  ordi_param block5;
} b6c96_params;

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
extern const float pass[2];
extern const float value[3];
extern const float scorevalue[6];
extern const float ownership[19][19];
extern const float ending[64];
extern const float input9[22][19][19];
extern const float policy9[2][19][19];
extern const float inputGlobal9[19];
extern const float afterblocks9[96][19][19];
extern const float after_policyrowsG9[96];
extern const float pass9[2];
extern const float value9[3];
extern const float scorevalue9[6];
extern const float ownership9[19][19];


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

float focus(float (*input)[3][3], int input_channel, float (*kernel)[3][3])
{
  float s = 0.0f;
  for(int i =0; i<input_channel; i++)
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

void fries(float (*input_padded)[21][21], int input_channel, float (*output)[3][3], int x, int y)
{
  for(int i=0; i<input_channel; i++)
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

void slip3x3(const float (*input)[19][19], int input_channel, float output_a_plane[19][19], float (*kernel)[3][3])
{
  // padding
  float (*input_padded)[19+1+1][19+1+1] = NULL;
  input_padded  = malloc(sizeof(float[19+1+1][19+1+1]) * input_channel);
  if (input_padded == NULL) { printf("malloc failed\n"); exit(1);}
  for(int i=0; i<input_channel; i++)
  {
    for(int j=0; j<19+1+1; j++)
    {
      for(int k=0; k<19+1+1; k++)
      {
        if(j==0 || k ==0 || j== 19+1 || k == 19+1)
          input_padded[i][j][k] = 0;
        else
          input_padded[i][j][k] = input[i][j-1][k-1];
      }
    }
  }
  // conv
  for(int i=0; i<19; i++)
  {
    for(int j=0; j<19; j++)
    {
      float (*fry)[3][3];
      fry = malloc(sizeof(float[3][3]) * input_channel);
      if (fry == NULL) { printf("malloc failed\n"); exit(1);}
      fries(input_padded, input_channel, fry, i, j);
      output_a_plane[i][j] = focus(fry, input_channel, kernel);
      free(fry);
    }
  }
  free(input_padded);
}

// for example: input[22][19][19], output[96][19][19], kernel[96][22][3][3]
void conv3x3(const float (*input)[19][19], int input_channel, float (*output)[19][19], int output_channel, float (*kernel)[3][3])
{
  for (int I=0; I < output_channel; I++)
  {
    float (*cannon)[3][3] = kernel + I * input_channel; /*for example: cannon[22][3][3]*/
    slip3x3(input, input_channel, output[I], cannon);
  }
}

void conv1x1(const float *input, int input_channel, float *output, int output_channel, int square_side, float *kernel)
{
  for (int I=0; I < output_channel; I++)
  {
    for(int i=0; i<square_side; i++)
      for(int j=0; j<square_side; j++)
      {
        float s = 0.0f;
        for (int k=0; k < input_channel; k++)
        {
          s+=input[k*square_side*square_side + i*square_side + j]*kernel[I*input_channel + k];
        }
        output[I*square_side*square_side + i*square_side + j] = s;
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
  
    if(diff > 0.001)
    {
      printf("err1 error: diff too big!! i=%d: output: %f, suppose: %f\n", i, mat1[i], mat2[i]);
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

void norm(const float (*input)[19][19], float (*output)[19][19], int channel, float *scale, float *bias, int board_size)
{
    for (int i = 0; i < channel; i++) {
      for (int j = 0; j < 19; ++j) {
        for (int k = 0; k < 19; ++k) {
          float x = input[i][j][k] * scale[i] + bias[i];
          if(j >= board_size || k >= board_size)
          {
            // out of board, mask filter
            output[i][j][k] = 0.0f;
          }
          else
          {
            output[i][j][k] = (x > 0.0f) ? x : 0.0f;
          }
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

void add_broadcast(const float (*input)[19][19], float (*output)[19][19], float *adder, int channel)
{
  for (int i = 0; i < channel; i++) {
    for (int j = 0; j < 19; ++j) {
      for (int k = 0; k < 19; ++k) {
        output[i][j][k] = input[i][j][k] + adder[i];
        }
      }
  }
}

void rowsG(float input[32][19][19], float output[96], int isValueHead, int board_size)
{
  float div = (float)(board_size * board_size);
  float sqrtdiv = (float)board_size;
  for(int I = 0; I<32; I++)
  {
    float s = 0.0f;
    float max = -1.0f;
    for(int i=0; i<19; i++)
      for(int j=0; j<19; j++)
      {
        float x = input[I][i][j];
        s += x;
        // katago原文
        // Init to -1.0 above and + mask - 1.0 is because it will effectively make all padded space into -1.0
        // which is lower than the lowest value that any current activation function will produce.
        // so the max over all valid spaces will the same as the mask over all spaces including padding
        // We're relying on all padded space being equal to 0 because this gpool only ever follows a BN+Activate with a mask.
        float maskVal = i<board_size && j<board_size? 1.0f : 0.0f;
        float temp = x + (maskVal - 1.0f);
        max = temp > max ? temp : max;
      }
        
    float mean = s / div;
    output[0 + I] = mean;
    output[32 + I] = mean * (sqrtdiv - 14.0f) * 0.1f;
    output[64 + I] = isValueHead? mean * ((sqrtdiv - 14.0f) * (sqrtdiv - 14.0f) * 0.01f - 0.1f) : max;
  }
}

void ordi(int board_size, float input[96][19][19], float output[96][19][19], float scale0[96], float bias0[96], const float (*afternorm)[19][19], float kernel1[96][96][3][3], const float (*afterconv)[19][19], float scale1[96], float bias1[96], const float (*afternorm2)[19][19], float kernel2[96][96][3][3], const float (*afterconv2)[19][19])
{
    float maxdiff;
    int erri, errj, errk;

    float output0[96][19][19];

    norm(input, output0, 96, scale0, bias0, board_size);

    if(afternorm != NULL)
    {
      printf("sumdiff: %f\n", err3((float*)output0, (const float*)afternorm, 96, 19, 19, &maxdiff, &erri, &errj, &errk));
      printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);
    }
    
    float output1[96][19][19];

    conv3x3(output0, 96, output1, 96, (float (*)[3][3])kernel1);

    if(afterconv != NULL)
    {
      printf("sumdiff: %f\n", err3((float*)output1, (const float*)afterconv, 96, 19, 19, &maxdiff, &erri, &errj, &errk));
      printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);
    }
    
    float output2[96][19][19];

    norm(output1, output2, 96, scale1, bias1, board_size);

    if(afternorm2 != NULL)
    {
      printf("sumdiff: %f\n", err3((float*)output2, (const float*)afternorm2, 96, 19, 19, &maxdiff, &erri, &errj, &errk));
      printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);
    }
    float output3[96][19][19];

    conv3x3(output2, 96, output3, 96, (float (*)[3][3])kernel2);

    add(input, output, output3);

    if(afterconv2 != NULL)
    {
      printf("sumdiff: %f\n", err3((float*)output, (const float*)afterconv2, 96, 19, 19, &maxdiff, &erri, &errj, &errk));
      printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);
    }

}

void gpool(int board_size, float input[96][19][19], float output[96][19][19], float scale0[96], float bias0[96], float kernel0[64][96][3][3], const float (*regularOut)[19][19], float kernel1[32][96][3][3], float scale1[32], float bias1[32], const float (*gpoolOut)[19][19], float nn[64][96], float scale2[64], float bias2[64], float kernel2[96][64][3][3], const float (*afterglobal)[19][19])
{
    float maxdiff;
    int erri, errj, errk;

    float output0[96][19][19];

    norm(input, output0, 96, scale0, bias0, board_size);

    //从这开始分为了两支，先是regular支

    float output1[64][19][19];

    conv3x3(output0, 96, output1, 64, (float (*)[3][3])kernel0);

    if(regularOut != NULL)
    {
      printf("sumdiff: %f\n", err3((float*)output1, (const float*)regularOut, 64, 19, 19, &maxdiff, &erri, &errj, &errk));
      printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);
    }

    // g分支

    float output2[32][19][19];

    conv3x3(output0, 96, output2, 32, (float (*)[3][3])kernel1);

    float output3[32][19][19];

    norm(output2, output3, 32, scale1, bias1, board_size);

    if(gpoolOut != NULL)
    {
      printf("sumdiff: %f\n", err3((float*)output3, (const float*)gpoolOut, 32, 19, 19, &maxdiff, &erri, &errj, &errk));
      printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);
    }
    
    float output4[96];
    rowsG(output3, output4, FALSE, board_size);

    float output5[64];
    conv1x1((float*)output4, 96, (float*)output5, 64, 1, (float*)nn);

    float output6[64][19][19];
    add_broadcast(output1, output6, output5, 64);

    // 汇聚后再来一次normconv
    // BINS[29]全是0，跳过
    // BINS[30]跳过

    float output7[64][19][19];

    norm(output6, output7, 64, scale2, bias2, board_size);

    float output8[96][19][19];

    conv3x3(output7, 64, output8, 96, (float (*)[3][3])kernel2);

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

void load_parameters_from_bin(b6c96_params* param)
{
  load_ordi(param->block0.scale0, param->block0.bias0, param->block0.kernel0, param->block0.scale1, param->block0.bias1, param->block0.kernel1, 3, 4, 5, 8, 9, 10);
  load_ordi(param->block1.scale0, param->block1.bias0, param->block1.kernel0, param->block1.scale1, param->block1.bias1, param->block1.kernel1, 12, 13, 14, 17, 18, 19);
  load_gpool(param->block2.scale0, param->block2.bias0, param->block2.kernel0, param->block2.kernel1, param->block2.scale1, param->block2.bias1, param->block2.nn0, param->block2.scale2, param->block2.bias2, param->block2.kernel2, 21, 22, 23, 24, 26, 27, 28, 31, 32, 33);
  load_ordi(param->block3.scale0, param->block3.bias0, param->block3.kernel0, param->block3.scale1, param->block3.bias1, param->block3.kernel1, 35, 36, 37, 40, 41, 42);
  load_gpool(param->block4.scale0, param->block4.bias0, param->block4.kernel0, param->block4.kernel1, param->block4.scale1, param->block4.bias1, param->block4.nn0, param->block4.scale2, param->block4.bias2, param->block4.kernel2, 44, 45, 46, 47, 49, 50, 51, 54, 55, 56);
  load_ordi(param->block5.scale0, param->block5.bias0, param->block5.kernel0, param->block5.scale1, param->block5.bias1, param->block5.kernel1, 58, 59, 60, 63, 64, 65);
}

void forward(int board_size, b6c96_params* param, const float input[22][19][19], const float inputGlobal[19], float policy[2][19][19], float pass[2], float value[3], float scorevalue[6], float ownership[19][19], const float (*afterblocks)[19][19], const float *after_policyrowsG9)
{
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
    conv3x3(input, 22, output0, 96, (float (*)[3][3])kernel0);

    float maxdiff;
    int erri, errj, errk;

    float output1[96];
    float linear0[96][19];
    n=0;

    for(int j=0; j <19; j++)
      for(int i=0; i < 96; i++)
      {
        linear0[i][j] = BINS[1].floats[n++];
      }

    conv1x1((float*)inputGlobal, 19, (float*)output1, 96, 1, (float*)linear0);

    float output2[96][19][19];
    add_broadcast(output0, output2, output1, 96);

    // block阶段
    // 6个block分别是 ordi, ordi, gpool, ordi, gpool, ordi

    float output3[96][19][19];
    float output4[96][19][19];
    float output5[96][19][19];
    float output6[96][19][19];
    float output7[96][19][19];
    float output8[96][19][19];

    ordi (board_size, output2, output3, param->block0.scale0,  param->block0.bias0,  NULL,     param->block0.kernel0,    NULL,      param->block0.scale1,  param->block0.bias1,  NULL,       param->block0.kernel1,  NULL);
    ordi (board_size, output3, output4, param->block1.scale0,  param->block1.bias0,  NULL,     param->block1.kernel0,    NULL,      param->block1.scale1,  param->block1.bias1,  NULL,       param->block1.kernel1,  NULL);
    gpool(board_size, output4, output5, param->block2.scale0,  param->block2.bias0,  param->block2.kernel0,   NULL,       param->block2.kernel1,   param->block2.scale1,  param->block2.bias1,  NULL,       param->block2.nn0,      param->block2.scale2,  param->block2.bias2,  param->block2.kernel2,  NULL);
    ordi (board_size, output5, output6, param->block3.scale0,  param->block3.bias0,  NULL,     param->block3.kernel0,    NULL,      param->block3.scale1,  param->block3.bias1,  NULL,       param->block3.kernel1,  NULL);
    gpool(board_size, output6, output7, param->block4.scale0,  param->block4.bias0,  param->block4.kernel0,   NULL,       param->block4.kernel1,   param->block4.scale1,  param->block4.bias1,  NULL,       param->block4.nn0,      param->block4.scale2,  param->block4.bias2,  param->block4.kernel2,  NULL);
    ordi (board_size, output7, output8, param->block5.scale0,  param->block5.bias0,  NULL,     param->block5.kernel0,    NULL,      param->block5.scale1,  param->block5.bias1,  NULL,       param->block5.kernel1,  NULL);

    if(afterblocks)
    {
      printf("sumdiff: %f\n", err3((float*)output8, (const float*)afterblocks, 96, 19, 19, &maxdiff, &erri, &errj, &errk));
      printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);
    }

    /* 6个block全部结束*/
    // 来一次final norm

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

    norm(output8, output9, 96, scale14, bias14, board_size);

    /* 以下开始各种头，policy, pass, value, scorevalue, ownership一共5个头 */

    /*分支1备用*/
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

    if(after_policyrowsG9)
    {
      printf("checking after_policyrowsG9 ...\n");
      printf("sumdiff: %f\n", err1((float*)output13, (const float*)after_policyrowsG9, 96, &maxdiff, &erri));
      printf("maxdiff: %f, %d\n", maxdiff, erri);
      printf("checking after_policyrowsG9 done.\n");
    }

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

    float nn3[32][96];
    n=0;

    for(int j=0; j < 96; j++)
      for(int i=0; i < 32; i++)
      {
        nn3[i][j] = BINS[79].floats[n];
        n++;
      }

    float output18[32];
    conv1x1((float*)output13, 96, (float*)output18, 32, 1, (float*)nn3);

    float adder0[32];
    n=0;
    for(int i=0; i < 32; i++)
    {
      adder0[i] = BINS[80].floats[n];
      n++;
    }
    float output19[32];
    // 加Bias
    for(int i=0; i < 32; i++)
    {
      output19[i] = output18[i] + adder0[i];
    }
    // 纯relu
    for(int i=0; i < 32; i++)
    {
      output19[i] = output19[i] > 0 ? output19[i] : 0.0f;
    }
    //norm
    float mult[2][32];
    n=0;

    for(int j=0; j < 32; j++)
      for(int i=0; i < 2; i++)
      {
        mult[i][j] = BINS[81].floats[n];
        n++;
      }
    
    for(int I=0; I < 2; I++)
    {
      float s = 0.0f;
      for(int i=0; i < 32; i++)
      {
        s += mult[I][i] * output19[i];
      }
      pass[I] = s;
    }

    /* Value 头*/
    //先Conv
    float kernel18[32][96];

    n=0;

    for(int j=0; j <96; j++)
      for(int i=0; i <32; i++)
      {
        kernel18[i][j] = BINS[82].floats[n++];
      }

    float output21[32][19][19];
    conv1x1((float*)output9, 96, (float*)output21, 32, 19, (float*)kernel18);

    //norm
    float scale17[32];
    float bias17[32];
    float output22[32][19][19];
    n=0;

    for(int i=0; i < 32; i++)
    {
      scale17[i] = BINS[84].floats[n];
      bias17[i] = BINS[85].floats[n];
      n++;
    }

    norm(output21, output22, 32, scale17, bias17, board_size);

    //ownership在此分出

    float output23[96];
    rowsG(output22, output23, TRUE, board_size);

    //mult
    float nn4[64][96];
    n=0;

    for(int j=0; j < 96; j++)
      for(int i=0; i < 64; i++)
      {
        nn4[i][j] = BINS[86].floats[n];
        n++;
      }

    float output24[64];
    conv1x1((float*)output23, 96, (float*)output24, 64, 1, (float*)nn4);

    float adder1[64];
    n=0;
    for(int i=0; i < 64; i++)
    {
      adder1[i] = BINS[87].floats[n];
      n++;
    }
    float output25[64];
    // 加Bias
    for(int i=0; i < 64; i++)
    {
      output25[i] = output24[i] + adder1[i];
    }
    // 纯relu
    for(int i=0; i < 64; i++)
    {
      output25[i] = output25[i] > 0 ? output25[i] : 0.0f;
    }

    //分出分支value和scorevalue

    /* value分支 */
    // mult
    float nn5[3][64];
    n=0;

    for(int j=0; j < 64; j++)
      for(int i=0; i < 3; i++)
      {
        nn5[i][j] = BINS[88].floats[n];
        n++;
      }

    float output26[3];
    conv1x1((float*)output25, 64, (float*)output26, 3, 1, (float*)nn5);

    float adder2[3];
    n=0;
    for(int i=0; i < 3; i++)
    {
      adder2[i] = BINS[89].floats[n];
      n++;
    }

    // 加Bias
    for(int i=0; i < 3; i++)
    {
      value[i] = output26[i] + adder2[i];
    }

    /* scorevalue分支 */
    float nn6[6][64];
    n=0;

    for(int j=0; j < 64; j++)
      for(int i=0; i < 6; i++)
      {
        nn6[i][j] = BINS[90].floats[n];
        n++;
      }
    
    float output28[6];
    conv1x1((float*)output25, 64, (float*)output28, 6, 1, (float*)nn6);

    float adder3[6];
    n=0;
    for(int i=0; i < 6; i++)
    {
      adder3[i] = BINS[91].floats[n];
      n++;
    }
    // 加Bias
    for(int i=0; i < 6; i++)
    {
      scorevalue[i] = output28[i] + adder3[i];
    }

    //conv

    // ownership
    float kernel19[32];

    n=0;

    for(int i=0; i <32; i++)
    {
      kernel19[i] = BINS[92].floats[n++];
    }

    conv1x1((float*)output22, 32, (float*)ownership, 1, 19, (float*)kernel19);

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

    printf("loading parameters...\n");
    b6c96_params *param = malloc(sizeof(b6c96_params));
    load_parameters_from_bin(param);
    printf("load parameters done.\n");


    printf("forwarding 19x19 ...\n");
    float out_policy[2][19][19], out_pass[2], out_value[3], out_scorevalue[6], out_ownership[19][19];
    forward(19, param, input, inputGlobal, out_policy, out_pass, out_value, out_scorevalue, out_ownership, afterblocks, NULL);
    printf("forwarding 19x19 done.\n");

    float maxdiff;
    int erri, errj, errk;

    printf("checking 19x19 policy ...\n");
    printf("sumdiff: %f\n", err3((float*)out_policy, (const float*)policy, 2, 19, 19, &maxdiff, &erri, &errj, &errk));
    printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);
    printf("check 19x19 policy done.\n");

    printf("checking 19x19 pass ...\n");
    printf("sumdiff: %f\n", err1((float*)out_pass, (const float*)pass, 2, &maxdiff, &erri));
    printf("maxdiff: %f, %d\n", maxdiff, erri);
    printf("check 19x19 pass done.\n");

    printf("checking 19x19 value ...\n");
    printf("sumdiff: %f\n", err1((float*)out_value, (const float*)value, 3, &maxdiff, &erri));
    printf("maxdiff: %f, %d\n", maxdiff, erri);
    printf("check 19x19 value done.\n");

    printf("checking 19x19 scorevalue ...\n");
    printf("sumdiff: %f\n", err1((float*)out_scorevalue, (const float*)scorevalue, 6, &maxdiff, &erri));
    printf("maxdiff: %f, %d\n", maxdiff, erri);
    printf("check 19x19 scorevalue done.\n");

    printf("checking 19x19 ownership ...\n");
    printf("sumdiff: %f\n", err3((float*)out_ownership, (const float*)ownership, 1, 19, 19, &maxdiff, &erri, &errj, &errk));
    printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);
    printf("check 19x19 ownership done.\n");

    printf("19x19 done.\n");


    /* 9x9 */

    printf("forwarding 9x9 ...\n");
    float out_policy9[2][19][19], out_pass9[2], out_value9[3], out_scorevalue9[6], out_ownership9[19][19];
    forward(9, param, input9, inputGlobal9, out_policy9, out_pass9, out_value9, out_scorevalue9, out_ownership9, afterblocks9, after_policyrowsG9);
    printf("forwarding 9x9 done.\n");

    for(int i=0;i<19;i++)
      for(int j=0;j<19;j++)
      {
        printf("%f, ", out_policy9[0][i][j]);
        if(j==18)
        {
          printf("\n");
        }
      }

    printf("checking 9x9 policy ...\n");
    printf("sumdiff: %f\n", err3((float*)out_policy9, (const float*)policy9, 2, 19, 19, &maxdiff, &erri, &errj, &errk));
    printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);
    printf("check 9x9 policy done.\n");

    printf("checking 9x9 pass ...\n");
    printf("sumdiff: %f\n", err1((float*)out_pass9, (const float*)pass9, 2, &maxdiff, &erri));
    printf("maxdiff: %f, %d\n", maxdiff, erri);
    printf("check 9x9 pass done.\n");

    printf("checking 9x9 value ...\n");
    printf("sumdiff: %f\n", err1((float*)out_value9, (const float*)value9, 3, &maxdiff, &erri));
    printf("maxdiff: %f, %d\n", maxdiff, erri);
    printf("check 9x9 value done.\n");

    printf("checking 9x9 scorevalue ...\n");
    printf("sumdiff: %f\n", err1((float*)out_scorevalue9, (const float*)scorevalue9, 6, &maxdiff, &erri));
    printf("maxdiff: %f, %d\n", maxdiff, erri);
    printf("check 9x9 scorevalue done.\n");

    printf("checking 9x9 ownership ...\n");
    printf("sumdiff: %f\n", err3((float*)out_ownership9, (const float*)ownership9, 1, 19, 19, &maxdiff, &erri, &errj, &errk));
    printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);
    printf("check 9x9 ownership done.\n");

    printf("19x19 done.\n");
    


    free(param);
    return 0;
}
