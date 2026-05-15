#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
 * Serial k-means on Iris dataset (k-meas/iris.txt)
 *
 * Data format (each line):
 *   sepal_len sepal_wid petal_len petal_wid label
 * We only use the first 4 features for clustering; last column is ignored.
 *
 * Clustering meaning (Iris):
 * - Features describe iris flower morphology (sepal/petal length/width, cm).
 * - k-means groups flowers by shape; with k=3 it often corresponds to the 3 species.
 *
 * Build & run:
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
    int n = read_iris(X);
    if (n <= 0) { fprintf(stderr, "No data loaded\n"); return 1; }

    /* init centroids: pick first K points (short & deterministic) */
    double C[K][DIM];
    for (int k = 0; k < K; k++) for (int j = 0; j < DIM; j++) C[k][j] = X[k][j];

    int assign[MAXN];
    for (int i = 0; i < n; i++) assign[i] = -1;

    const int maxIter = 100;
    for (int it = 0; it < maxIter; it++) {
        int changed = 0;

        /* assignment step */
        for (int i = 0; i < n; i++) {
            int bestk = 0;
            double bestd = dist2(X[i], C[0]);
            for (int k = 1; k < K; k++) {
                double d = dist2(X[i], C[k]);
                if (d < bestd) { bestd = d; bestk = k; }
            }
            if (assign[i] != bestk) { assign[i] = bestk; changed = 1; }
        }

        /* update step */
        double sum[K][DIM]; int cnt[K];
        memset(sum, 0, sizeof(sum));
        memset(cnt, 0, sizeof(cnt));
        for (int i = 0; i < n; i++) {
            int k = assign[i];
            cnt[k]++;
            for (int j = 0; j < DIM; j++) sum[k][j] += X[i][j];
        }
        for (int k = 0; k < K; k++) {
            if (cnt[k] == 0) continue; /* empty cluster: keep centroid */
            for (int j = 0; j < DIM; j++) C[k][j] = sum[k][j] / cnt[k];
        }

        if (!changed) break;
    }

    /* report */
    int cnt[K]; memset(cnt, 0, sizeof(cnt));
    double sse = 0.0;
    for (int i = 0; i < n; i++) { cnt[assign[i]]++; sse += dist2(X[i], C[assign[i]]); }

    printf("Serial k-means on iris: n=%d, k=%d, dim=%d\n", n, K, DIM);
    for (int k = 0; k < K; k++) {
        printf("C%d (count=%d):", k, cnt[k]);
        for (int j = 0; j < DIM; j++) printf(" %.4f", C[k][j]);
        printf("\n");
    }
    printf("SSE=%.6f\n", sse);
    return 0;
}
