#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "vector_matrix.h"


double cal_error(int m, double *beta, double *x, double *y);

void order(int m, int n, double **A, double *x, double *y);

void cal_centroid(int n, double **A, double *c);

void construct_A(int n, double *beta, double **A) {
    // 随机生成测试点，构造测试点矩阵A
    int i, j;
    double range = 0.05;

    for (i = 0; i < n + 1; i++) {
        for (j = 0; j < n; j++) {
            A[i][j] = ((rand() % ((int) (2 * range * 1e6 + 1))) + (beta[j] - range) * 1e6) / 1e6;
        }
    }
}

void reflection(int n, double **A, double *xo, double *xr) {
    double alpha = 1.0;
    double *t;
    int i;

    t = (double *) malloc(sizeof(double) * n);

    for (i = 0; i < n; i++) {
        xr[i] = 0;
    }
    vector_minus_vector(n, xo, A[n], t);
    num_times_vector(n, t, alpha, t);
    vector_add_vector(n, xo, t, xr);

    free(t);
}

void expansion(int n, double *xo, double *xr, double *xe) {
    double gama = 2.0;
    double *t;
    int i;

    t = (double *) malloc(sizeof(double) * n);

    for (i = 0; i < n; i++) {
        xe[i] = 0;
    }
    vector_minus_vector(n, xr, xo, t);
    num_times_vector(n, t, gama, t);
    vector_add_vector(n, xo, t, xe);

    free(t);
}

void shrink(int n, double **A) {
    double sigma = 0.5;
    double *t;
    int i;

    t = (double *) malloc(sizeof(double) * n);

    for (i = 1; i < n + 1; i++) {
        vector_minus_vector(n, A[i], A[0], t);
        num_times_vector(n, t, sigma, t);
        vector_add_vector(n, A[0], t, A[i]);
    }

    free(t);
}

void contraction(int m, int n, double **A, double *xo, double *xr, double *xc, double *x, double *y) {
    double rho = 0.5;
    double *t;
    int i;

    t = (double *) malloc(sizeof(double) * n);

    for (i = 0; i < n; i++) {
        xc[i] = 0;
    }
    if (cal_error(m, xr, x, y) < cal_error(m, A[n], x, y)) {
        printf("f(xr)<f(xn+1)\n");
        vector_minus_vector(n, xr, xo, t);
        num_times_vector(n, t, rho, t);
        vector_add_vector(n, xo, t, xc);
        if (cal_error(m, xc, x, y) < cal_error(m, xr, x, y)) {
            copy_vector(n, xc, A[n]);
            printf("f(xc)<f(xr)\n");
        } else {
            shrink(n, A);
            printf("f(xc)<f(xr)\n");
        }
    } else {
        printf("f(xr)>=f(xn+1)\n");
        vector_minus_vector(n, A[n], xo, t);
        num_times_vector(n, t, rho, t);
        vector_add_vector(n, xo, t, xc);
        if (cal_error(m, xc, x, y) < cal_error(m, A[n], x, y)) {
            copy_vector(n, xc, A[n]);
            printf("f(xc)<f(xn+1)\n");
        } else {
            shrink(n, A);
            printf("f(xc)>=f(xn+1)\n");
        }
    }

    free(t);
}

void NMMethod(int m, int n, double **A, double *beta, double *x, double *y, double eps) {
    double *xo, *xr, *xe, *xc;
    double e = 1.0;
    int i, j, step = 0;

    xo = (double *) malloc(sizeof(double) * n);
    xr = (double *) malloc(sizeof(double) * n);
    xe = (double *) malloc(sizeof(double) * n);
    xc = (double *) malloc(sizeof(double) * n);

    construct_A(n, beta, A);
    printf("A = \n");
    for (i = 0; i < n + 1; i++) {
        for (j = 0; j < n; j++) {
            printf("%lf ", A[i][j]);
        }
        printf("\n");
    }
    order(m, n, A, x, y);

    while (e > eps && fabs((cal_error(m, A[0], x, y) - cal_error(m, A[n], x, y))) > 1e-14) {
        printf("step = %d\n", step);
        step += 1;
        order(m, n, A, x, y);

        printf("A[0] = \n");
        for (i = 0; i < n; i++) {
            printf("%.10lf ", A[0][i]);
        }

        e = cal_error(m, A[0], x, y);
        printf("\ne = %lf\n", e);
        if (e < eps) {
            break;
        }

        cal_centroid(n, A, xo);
        reflection(n, A, xo, xr);

        if (cal_error(m, A[0], x, y) <= cal_error(m, xr, x, y) &&
            cal_error(m, xr, x, y) < cal_error(m, A[n - 1], x, y)) {
            copy_vector(n, xr, A[n]);
            printf("f(x1)<=f(xr)<f(xn)\n");
            continue;
        } else if (cal_error(m, xr, x, y) < cal_error(m, A[0], x, y)) {
            printf("f(xr)<f(x1)\n");
            expansion(n, xo, xr, xe);
            if (cal_error(m, xe, x, y) < cal_error(m, xr, x, y)) {
                copy_vector(n, xe, A[n]);
                printf("f(xe)<f(xr)\n");
                continue;
            } else {
                printf("f(xe)>=f(xr)\n");
                copy_vector(n, xr, A[n]);
                continue;
            }
        } else {
            printf("f(xr)>=f(xn)\n");
            contraction(m, n, A, xo, xr, xc, x, y);
        }
    }

    free(xo);
    free(xr);
    free(xe);
    free(xc);
}

int main() {
    /*
     该程序用于使用Nelder–Mead算法解决非线性最小二乘问题
     p(x)=a_0+a_1*sin(a_2*x)+a_3*cos(a_4*x)

     m 表示数据数目
     n 表示参数数目
     x 表示拟合点坐标
     y 表示拟合点函数值
     beta 表示参数向量
     A 为(n+1)*n阶测试点矩阵
     eps 表示可允许误差
     e 表示误差
    */
    int m = 11, n = 5;
    double x[] = {-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    double y[] = {-0.8669, -0.2997, 0.3430, 1.0072, 1.6416, 2.2022, 2.6558, 2.9823, 3.1755, 3.2416, 3.1974};
    // double beta[] = {1.3, 2.2, 1.1, 0.9, 1.8};
    double beta[] = {0.1, 0.1, 0.1, 0.1, 0.1};
    double **A;
    double eps = 0.00001;
    int i, j;
    srand((unsigned) time(NULL));

    A = (double **) malloc(sizeof(double *) * (n + 1));
    for (i = 0; i < n + 1; i++) {
        A[i] = (double *) malloc(sizeof(double) * n);
    }

    NMMethod(m, n, A, beta, x, y, eps);

    free(A);
}


double cal_error(int m, double *beta, double *x, double *y) {
    // 计算误差
    int i, j;
    double e = 0.0, t;

    for (i = 0; i < m; i++) {
        t = beta[0] + beta[1] * sin(beta[2] * x[i]) + beta[3] * cos(beta[4] * x[i]);
        e += pow(t - y[i], 2);
    }

    return sqrt(e);
}

void order(int m, int n, double **A, double *x, double *y) {
    // 对A中向量进行排序
    int i, j, k;
    double *t;

    t = (double *) malloc(n * sizeof(double));
    for (i = 0; i < n; i++) {
        k = i;
        for (j = i + 1; j < n + 1; j++) {
            if (cal_error(m, A[j], x, y) <= cal_error(m, A[k], x, y)) {
                k = j;
            }
        }
        if (k != i) {
            copy_vector(n, A[i], t);
            copy_vector(n, A[k], A[i]);
            copy_vector(n, t, A[k]);
        }
    }
    free(t);
}

void cal_centroid(int n, double **A, double *c) {
    // 计算单纯形几何中心坐标
    int i;

    for (i = 0; i < n; i++) {
        c[i] = 0;
    }
    for (i = 0; i < n; i++) {
        vector_add_vector(n, A[i], c, c);
    }
    vector_divide_num(n, c, n, c);
}