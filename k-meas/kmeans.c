#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
 * 在 Iris 数据集上的串行 k-means 算法 (k-meas/iris.txt)
 *
 * 数据格式（每行）：
 *   花萼长度 花萼宽度 花瓣长度 花瓣宽度 标签
 * 我们只使用前 4 个特征进行聚类；最后一列标签被忽略。
 *
 * 聚类意义 (Iris)：
 * - 特征描述了鸢尾花的形态（花萼/花瓣的长度/宽度，单位 cm）。
 * - k-means 根据形状对花进行分组；当 k=3 时，通常对应 3 个品种。
 *
 * 编译与运行：
 *   gcc -O2 -o kmeans k-meas/kmeans.c -lm
 *   ./kmeans
 */

#define DIM 4
#define K 3
#define MAXN 200

static int read_iris(double X[MAXN][DIM]) {
    FILE *fp = fopen("k-meas/iris.txt", "r");
    if (!fp) { perror("open iris.txt"); exit(1); }

    int n = 0;
    char line[256];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue;
        if (n >= MAXN) break;

        double a, b, c, d;
        int label;
        if (sscanf(line, "%lf %lf %lf %lf %d", &a, &b, &c, &d, &label) == 5) {
            X[n][0] = a; X[n][1] = b; X[n][2] = c; X[n][3] = d;
            n++;
        }
    }
    fclose(fp);
    return n;
}

static double dist2(const double *x, const double *c) {
    double s = 0.0;
    for (int j = 0; j < DIM; j++) { double t = x[j] - c[j]; s += t * t; }
    return s;
}

int main(void) {
    double X[MAXN][DIM];
    int n = read_iris(X); // 读取Iris数据集
    if (n <= 0) { fprintf(stderr, "No data loaded\n"); return 1; }

    /* 初始化聚类中心：选取前 K 个点（简单且确定） */
    double C[K][DIM];
    for (int k = 0; k < K; k++) for (int j = 0; j < DIM; j++) C[k][j] = X[k][j]; // 初始聚类中心定为前K个点

    int assign[MAXN];
    for (int i = 0; i < n; i++) assign[i] = -1;

    const int maxIter = 100;
    for (int it = 0; it < maxIter; it++) {
        int changed = 0;

        /* 分配步骤 */
        for (int i = 0; i < n; i++) {
            int bestk = 0;
            double bestd = dist2(X[i], C[0]);
            for (int k = 1; k < K; k++) {
                double d = dist2(X[i], C[k]);
                if (d < bestd) { bestd = d; bestk = k; } // 寻找最近的聚类中心
            }
            if (assign[i] != bestk) { assign[i] = bestk; changed = 1; }
        }

        /* 更新步骤 */
        double sum[K][DIM]; int cnt[K];
        memset(sum, 0, sizeof(sum));
        memset(cnt, 0, sizeof(cnt));
        for (int i = 0; i < n; i++) {
            int k = assign[i];
            cnt[k]++; // 统计属于该簇的点数
            for (int j = 0; j < DIM; j++) sum[k][j] += X[i][j]; // 累加坐标以便计算均值
        }
        for (int k = 0; k < K; k++) {
            if (cnt[k] == 0) continue;
            for (int j = 0; j < DIM; j++) C[k][j] = sum[k][j] / cnt[k]; // 更新聚类中心为均值点
        }

        if (!changed) break; // 如果所有点分配未变化，说明已收敛
    }

    /* 报告结果 */
    int cnt[K]; memset(cnt, 0, sizeof(cnt));
    double sse = 0.0;
    for (int i = 0; i < n; i++) { cnt[assign[i]]++; sse += dist2(X[i], C[assign[i]]); } // 计算误差平方和SSE

    printf("Serial k-means on iris: n=%d, k=%d, dim=%d\n", n, K, DIM);
    for (int k = 0; k < K; k++) {
        printf("C%d (count=%d):", k, cnt[k]);
        for (int j = 0; j < DIM; j++) printf(" %.4f", C[k][j]);
        printf("\n");
    }
    printf("SSE=%.6f\n", sse);
    return 0;
}
