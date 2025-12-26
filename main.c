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
extern const float after_second_norm_norm[96][19][19];
extern const float after2norms[96][19][19];
extern const float after2ordis[96][19][19];
extern const float regularOut[64][19][19];
extern const float gpoolOut2[32][19][19];

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
    
    // BINS[2]全是0，也没用到，跳过


    float output3[96][19][19];
    float scale0[96];
    float bias0[96];
    n=0;

    for(int i=0; i < 96; i++)
    {
      scale0[i] = BINS[3].floats[n];
      bias0[i] = BINS[4].floats[n];
      n++;
    }

    norm(output2, output3, scale0, bias0);
    printf("sumdiff: %f\n", err3((float*)output3, (const float*)afternorm, 96, 19, 19, &maxdiff, &erri, &errj, &errk));

    printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);


    float output4[96][19][19];
    float kernel1[96][96][3][3];
    n=0;

    for(int k=0; k <3; k++)
      for(int l=0; l <3; l++)
        for(int j=0; j <96; j++)
          for(int i=0; i <96; i++)
          {
            kernel1[i][j][k][l] = BINS[5].floats[n++];
          }

    conv9696(output3, output4, kernel1);
    printf("sumdiff: %f\n", err3((float*)output4, (const float*)afterconv, 96, 19, 19, &maxdiff, &erri, &errj, &errk));

    printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);


    // BINS[6]全是0，也没用到，跳过
    // BINS[7]全是1，但此时scale不像上次一样，没用这个。而是跳过了这个，用了下一个。

    float output5[96][19][19];
    float scale1[96];
    float bias1[96];
    n=0;

    for(int i=0; i < 96; i++)
    {
      scale1[i] = BINS[8].floats[n];
      bias1[i] = BINS[9].floats[n];
      n++;
    }

    norm(output4, output5, scale1, bias1);

    printf("sumdiff: %f\n", err3((float*)output5, (const float*)after_second_norm_norm, 96, 19, 19, &maxdiff, &erri, &errj, &errk));

    printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);

    float output6[96][19][19];
    float kernel2[96][96][3][3];
    n=0;

    for(int k=0; k <3; k++)
      for(int l=0; l <3; l++)
        for(int j=0; j <96; j++)
          for(int i=0; i <96; i++)
          {
            kernel2[i][j][k][l] = BINS[10].floats[n++];
          }

    conv9696(output5, output6, kernel2);

    float output7[96][19][19];

    add(output2, output7, output6);

    printf("sumdiff: %f\n", err3((float*)output7, (const float*)after2norms, 96, 19, 19, &maxdiff, &erri, &errj, &errk));

    printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);

    /*第二个block开始，还是ordi*/

    float output8[96][19][19];
    float scale2[96];
    float bias2[96];
    n=0;

    for(int i=0; i < 96; i++)
    {
      scale2[i] = BINS[12].floats[n];
      bias2[i] = BINS[13].floats[n];
      n++;
    }

    norm(output7, output8, scale2, bias2);

    float output9[96][19][19];
    float kernel3[96][96][3][3];
    n=0;

    for(int k=0; k <3; k++)
      for(int l=0; l <3; l++)
        for(int j=0; j <96; j++)
          for(int i=0; i <96; i++)
          {
            kernel3[i][j][k][l] = BINS[14].floats[n++];
          }

    conv9696(output8, output9, kernel3);

    // BINS[16]全是0，也没用到，跳过
    // BINS[17]全是1，但此时scale不像上次一样，没用这个。而是跳过了这个，用了下一个。

    float output10[96][19][19];
    float scale3[96];
    float bias3[96];
    n=0;

    for(int i=0; i < 96; i++)
    {
      scale3[i] = BINS[17].floats[n];
      bias3[i] = BINS[18].floats[n];
      n++;
    }

    norm(output9, output10, scale3, bias3);

    float output11[96][19][19];
    float kernel4[96][96][3][3];
    n=0;

    for(int k=0; k <3; k++)
      for(int l=0; l <3; l++)
        for(int j=0; j <96; j++)
          for(int i=0; i <96; i++)
          {
            kernel4[i][j][k][l] = BINS[19].floats[n++];
          }

    conv9696(output10, output11, kernel4);

    float output12[96][19][19];

    add(output7, output12, output11);

    printf("sumdiff: %f\n", err3((float*)output12, (const float*)after2ordis, 96, 19, 19, &maxdiff, &erri, &errj, &errk));

    printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);

    // 到了第三个block，gpool

    // BINS[20]全是0，跳过

    float output13[96][19][19];
    float scale4[96];
    float bias4[96];
    n=0;

    for(int i=0; i < 96; i++)
    {
      scale4[i] = BINS[21].floats[n];
      bias4[i] = BINS[22].floats[n];
      n++;
    }

    norm(output12, output13, scale4, bias4);

    //从这开始分为了两支，先是regular支

    float output14[64][19][19];
    float kernel5[64][96][3][3];
    n=0;

    for(int k=0; k <3; k++)
      for(int l=0; l <3; l++)
        for(int j=0; j <96; j++)
          for(int i=0; i <64; i++)
          {
            kernel5[i][j][k][l] = BINS[23].floats[n++];
          }

    conv9664(output13, output14, kernel5);

    printf("sumdiff: %f\n", err3((float*)output14, (const float*)regularOut, 64, 19, 19, &maxdiff, &erri, &errj, &errk));

    printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);

    // g分支

    float output15[32][19][19];
    float kernel6[32][96][3][3];
    n=0;

    for(int k=0; k <3; k++)
      for(int l=0; l <3; l++)
        for(int j=0; j <96; j++)
          for(int i=0; i <32; i++)
          {
            kernel6[i][j][k][l] = BINS[24].floats[n++];
          }

    conv9632(output13, output15, kernel6);

    float output16[32][19][19];
    float scale5[32];
    float bias5[32];
    n=0;

    for(int i=0; i < 32; i++)
    {
      scale5[i] = BINS[26].floats[n];
      bias5[i] = BINS[27].floats[n];
      n++;
    }

    norm32(output15, output16, scale5, bias5);

    printf("sumdiff: %f\n", err3((float*)output16, (const float*)gpoolOut2, 32, 19, 19, &maxdiff, &erri, &errj, &errk));

    printf("maxdiff: %f, %d, %d, %d\n", maxdiff, erri, errj, errk);










    return 0;
}
