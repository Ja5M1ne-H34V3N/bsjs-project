#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
 * 在 Iris 数据集上的 MPI 并行 k-means 算法 (k-meas/iris.txt)
 *
 * 数据集意义：
 * - 每个样本是一朵鸢尾花，由 4 个形态特征（单位 cm）描述：
 *   花萼长度、花萼宽度、花瓣长度、花瓣宽度。
 * - 聚类根据形状对花进行分组；当 k=3 时，通常与 3 个品种对齐。
 *
 * 并行思路 (重点说明)：
 * - 不是“每个聚类中心一个进程”（K 很小，这样会浪费并行性）。
 * - 我们采用“数据并行”：
 *   每个 MPI 进程分配约等量的数据点，计算这些点到所有 K 个聚类中心的距离，
 *   累加每个簇的局部总和/计数；然后通过 MPI_Allreduce 汇总形成新的聚类中心。
 *
 * 编译与运行：
 *   mpicc -O2 -o mpikmeans k-meas/mpikmeans.c -lm
 *   mpirun -np 4 ./mpikmeans
 */

#define DIM 4
#define K 3
#define MAXN 200

static int read_iris(double X[MAXN][DIM], int rank) {
    FILE *fp = fopen("k-meas/iris.txt", "r");
    if (!fp) { if(rank==0) perror("open iris.txt"); return -1; }

    int n = 0;
    char line[256];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue;
        if (n >= MAXN) break;
        double a,b,c,d; int label;
        if (sscanf(line, "%lf %lf %lf %lf %d", &a,&b,&c,&d,&label) == 5) {
            X[n][0]=a; X[n][1]=b; X[n][2]=c; X[n][3]=d;
            n++;
        }
    }
    fclose(fp);
    return n;
}

static double dist2(const double *x, const double *c) {
    double s=0.0;
    for (int j=0;j<DIM;j++){ double t=x[j]-c[j]; s+=t*t; }
    return s;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    /* 根进程读取数据，然后广播所有数据（对于小规模的 Iris 数据集来说简单快捷） */
    double X[MAXN][DIM];
    int n = 0;
    if (rank == 0) {
        n = read_iris(X, rank);
        if (n <= 0) { fprintf(stderr,"No data loaded\n"); MPI_Abort(MPI_COMM_WORLD,1); }
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); // 广播数据总量n到所有进程
    MPI_Bcast(X, MAXN*DIM, MPI_DOUBLE, 0, MPI_COMM_WORLD); // 广播整个数据集X

    /* 数据划分: [begin, end) */
    int begin = (n * rank) / size; // 计算当前进程处理数据的起始索引
    int end   = (n * (rank + 1)) / size; // 计算当前进程处理数据的结束索引

    /* 初始化聚类中心：使用前 K 个点（确定性） */
    double C[K][DIM];
    if (rank == 0) {
        for (int k=0;k<K;k++) for (int j=0;j<DIM;j++) C[k][j] = X[k][j];
    }
    MPI_Bcast(C, K*DIM, MPI_DOUBLE, 0, MPI_COMM_WORLD); // 广播初始聚类中心

    int assign[MAXN];
    for (int i=begin;i<end;i++) assign[i] = -1;

    const int maxIter = 100;
    for (int it=0; it<maxIter; it++) {
        int changed_local = 0;

        /* 本地分配 + 累加本地总和/计数 */
        double sum_local[K][DIM]; int cnt_local[K];
        memset(sum_local, 0, sizeof(sum_local));
        memset(cnt_local, 0, sizeof(cnt_local));

        for (int i=begin;i<end;i++) {
            int bestk = 0;
            double bestd = dist2(X[i], C[0]);
            for (int k=1;k<K;k++) {
                double d = dist2(X[i], C[k]);
                if (d < bestd) { bestd=d; bestk=k; }
            }
            if (assign[i] != bestk) { assign[i]=bestk; changed_local = 1; } // 标记本地分配是否有变化

            cnt_local[bestk]++; // 累计本地属于该簇的点数
            for (int j=0;j<DIM;j++) sum_local[bestk][j] += X[i][j]; // 累加本地部分和
        }

        /* 全局规约以获取新的聚类中心 */
        double sum_global[K][DIM]; int cnt_global[K];
        MPI_Allreduce(sum_local, sum_global, K*DIM, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // 规约所有进程的部分和
        MPI_Allreduce(cnt_local, cnt_global, K,     MPI_INT,    MPI_SUM, MPI_COMM_WORLD); // 规约所有进程的计数

        for (int k=0;k<K;k++) {
            if (cnt_global[k] == 0) continue;
            for (int j=0;j<DIM;j++) C[k][j] = sum_global[k][j] / cnt_global[k]; // 更新全局聚类中心
        }

        int changed_global = 0;
        MPI_Allreduce(&changed_local, &changed_global, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD); // 汇集所有进程的变更状态，只要有一个进程变了就继续迭代
        if (!changed_global) break;
    }

    /* 报告结果 (根进程): 簇大小和 SSE */
    int cnt_local[K]; memset(cnt_local,0,sizeof(cnt_local));
    double sse_local = 0.0;
    for (int i=begin;i<end;i++){
        cnt_local[assign[i]]++;
        sse_local += dist2(X[i], C[assign[i]]);
    }
    int cnt_global[K]; double sse_global = 0.0;
    MPI_Reduce(cnt_local, cnt_global, K, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); // 汇总最终计数到主进程
    MPI_Reduce(&sse_local, &sse_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 汇总最终SSE到主进程

    if (rank == 0) {
        printf("MPI k-means on iris: n=%d, k=%d, dim=%d, P=%d\n", n, K, DIM, size);
        for (int k=0;k<K;k++) {
            printf("C%d (count=%d):", k, cnt_global[k]);
            for (int j=0;j<DIM;j++) printf(" %.4f", C[k][j]);
            printf("\n");
        }
        printf("SSE=%.6f\n", sse_global);
        printf("Parallel strategy: data-parallel (each rank handles ~n/P points, Allreduce sums/counts)\n");
    }

    MPI_Finalize();
    return 0;
}
