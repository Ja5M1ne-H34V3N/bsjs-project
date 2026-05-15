#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
 * MPI parallel k-means on Iris dataset (k-meas/iris.txt)
 *
 * Dataset meaning:
 * - Each sample is an iris flower described by 4 morphology features (cm):
 *   sepal_length, sepal_width, petal_length, petal_width.
 * - Clustering groups flowers by shape; with k=3 it often aligns with 3 species.
 *
 * Parallel idea (重点说明):
 * - NOT "one process per centroid" (K is small, would waste parallelism).
 * - We do "data parallelism":
 *   each MPI rank gets ~equal number of data points, computes distances to all K centroids,
 *   accumulates partial sums/counts for each cluster; then MPI_Allreduce to form new centroids.
 *
 * Build & run:
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

    /* root reads, then broadcast all data (simple & short for small Iris dataset) */
    double X[MAXN][DIM];
    int n = 0;
    if (rank == 0) {
        n = read_iris(X, rank);
        if (n <= 0) { fprintf(stderr,"No data loaded\n"); MPI_Abort(MPI_COMM_WORLD,1); }
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(X, MAXN*DIM, MPI_DOUBLE, 0, MPI_COMM_WORLD); /* n<=MAXN; broadcast full buffer for simplicity */

    /* data partition: [begin, end) */
    int begin = (n * rank) / size;
    int end   = (n * (rank + 1)) / size;

    /* init centroids: use first K points (deterministic) */
    double C[K][DIM];
    if (rank == 0) {
        for (int k=0;k<K;k++) for (int j=0;j<DIM;j++) C[k][j] = X[k][j];
    }
    MPI_Bcast(C, K*DIM, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int assign[MAXN];
    for (int i=begin;i<end;i++) assign[i] = -1;

    const int maxIter = 100;
    for (int it=0; it<maxIter; it++) {
        int changed_local = 0;

        /* local assignment + accumulate local sums/counts */
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
            if (assign[i] != bestk) { assign[i]=bestk; changed_local = 1; }

            cnt_local[bestk]++;
            for (int j=0;j<DIM;j++) sum_local[bestk][j] += X[i][j];
        }

        /* global reduce to get new centroids */
        double sum_global[K][DIM]; int cnt_global[K];
        MPI_Allreduce(sum_local, sum_global, K*DIM, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(cnt_local, cnt_global, K,     MPI_INT,    MPI_SUM, MPI_COMM_WORLD);

        for (int k=0;k<K;k++) {
            if (cnt_global[k] == 0) continue;
            for (int j=0;j<DIM;j++) C[k][j] = sum_global[k][j] / cnt_global[k];
        }

        int changed_global = 0;
        MPI_Allreduce(&changed_local, &changed_global, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        if (!changed_global) break;
    }

    /* report (root): cluster sizes and SSE */
    int cnt_local[K]; memset(cnt_local,0,sizeof(cnt_local));
    double sse_local = 0.0;
    for (int i=begin;i<end;i++){
        cnt_local[assign[i]]++;
        sse_local += dist2(X[i], C[assign[i]]);
    }
    int cnt_global[K]; double sse_global = 0.0;
    MPI_Reduce(cnt_local, cnt_global, K, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sse_local, &sse_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

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
