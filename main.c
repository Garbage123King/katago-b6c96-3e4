#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

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

    printf("解压完成：%s\n", binfile);

    /* 4. 打开解压出的 .bin 文件，按 float32 读取 */
    FILE *fp = fopen(binfile, "rb");
    if (!fp) {
        fprintf(stderr, "无法打开 %s\n", binfile);
        return 1;
    }

    float value;
    size_t count = 0;

    printf("读取浮点数：\n");

    while (fread(&value, sizeof(float), 1, fp) == 1) {
        printf("%zu: %f\n", count, value);
        count++;
    }

    fclose(fp);

    printf("总共读取 %zu 个 float32\n", count);

    return 0;
}
