#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector_matrix.h"

void trace(int n, double **A, double **trA);

void output(int m, int n, double *beta);

void construct_J(int m, double *beta, double *x, double **J) {
    // 构造Jacobi矩阵
    int i;

    for (i = 0; i < m; i++) {
        J[i][0] = 1.0;
    }
    for (i = 0; i < m; i++) {
        J[i][1] = sin(beta[2] * x[i]);
    }
    for (i = 0; i < m; i++) {
        J[i][2] = beta[1] * x[i] * cos(beta[2] * x[i]);
    }
    for (i = 0; i < m; i++) {
        J[i][3] = cos(beta[4] * x[i]);
    }
    for (i = 0; i < m; i++) {
        J[i][4] = (-beta[3] * x[i] * sin(beta[4] * x[i]));
    }
}

void construct_r(int m, double *x, double *y, double *beta, double *r) {
    // 构造r
    int i;

    for (i = 0; i < m; i++) {
        r[i] = y[i] - (beta[0] + beta[1] * sin(beta[2] * x[i]) + beta[3] * cos(beta[4] * x[i]));
    }
}


void LMAlgorithm(int m, int n, double *beta, double alpha, double **J, double *r, double lambda) {
    /*
     JT 为J的转置
     JTJ 为JT*J
     JTr 为JT*r
     delta_beta 为Δ
     trJTJ 为对角线元素为JTJ对角线元素，其他元素均为0的矩阵
     */
    double **JT, **JTJ, **trJTJ;
    double *JTr, *delta_beta;
    int i, j;

    JT = (double **) malloc(sizeof(double *) * n);
    for (i = 0; i < n; i++) {
        JT[i] = (double *) malloc(sizeof(double) * m);
    }
    JTJ = (double **) malloc(sizeof(double *) * n);
    for (i = 0; i < n; i++) {
        JTJ[i] = (double *) malloc(sizeof(double) * n);
    }
    trJTJ = (double **) malloc(sizeof(double *) * n);
    for (i = 0; i < n; i++) {
        trJTJ[i] = (double *) malloc(sizeof(double) * n);
    }
    JTr = (double *) malloc(n * sizeof(double));
    delta_beta = (double *) malloc(n * sizeof(double));

    matrix_tranform(m, n, J, JT);
    matrix_times_matrix(n, m, n, JT, J, JTJ);
    trace(n, JTJ, trJTJ);
    num_times_matrix(n, n, trJTJ, lambda, trJTJ);
    matrix_add_matrix(n, n, JTJ, trJTJ, JTJ);

    matrix_times_vector(n, m, JT, r, JTr);

    solve_linear(n, JTJ, JTr, delta_beta);
    num_times_vector(n, delta_beta, alpha, delta_beta);
    vector_add_vector(n, beta, delta_beta, beta);

    free(JT);
    free(JTJ);
    free(trJTJ);
    free(JTr);
    free(delta_beta);
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


int main() {
    /*
     该程序用于使用Levenberg–Marquardt算法解决非线性最小二乘问题
     p(x)=a_0+a_1*sin(a_2*x)+a_3*cos(a_4*x)

     m 表示数据数目
     n 表示参数数目
     x 表示拟合点坐标
     y 表示拟合点函数值
     beta 表示参数向量
     J 为Jacobi矩阵
     r 为残差
     alpha 表示梯度下降速率
     lambda 表示阻尼因子
     nu1,nu2 表示迭代因子
     eps 表示可允许误差
     e 表示误差
     d 表示前后两次误差之差
     delta 表示d可允许最大值
    */
    int m = 11, n = 5;
    double x[] = {-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    double y[] = {-0.8669, -0.2997, 0.3430, 1.0072, 1.6416, 2.2022, 2.6558, 2.9823, 3.1755, 3.2416, 3.1974};
    double beta[] = {0.1, 0.1, 0.1, 0.1, 0.1};
    // double beta[] = {1.3, 2.2, 1.1, 0.9, 1.8};
    double **J;
    double *r;
    double alpha = 0.001, eps = 0.00001, e = 10, delta = 1e-16, d = 1, e_ = 0;
    double lambda = 0.1, nu1 = 2.0, nu2 = 3.0;
    int step = 0;
    int i, j;

    // J是一个m*n矩阵
    J = (double **) malloc(sizeof(double *) * m);
    for (i = 0; i < m; i++) {
        J[i] = (double *) malloc(sizeof(double) * n);
    }
    r = (double *) malloc(m * sizeof(double));

    while (e > eps && fabs(d) > delta) {
        construct_J(m, beta, x, J);
        construct_r(m, x, y, beta, r);
        LMAlgorithm(m, n, beta, alpha, J, r, lambda);
        e = cal_error(m, beta, x, y);
        output(m, n, beta);
        d = e_ - e;
        e_ = e;
        if (d > 0) {
            // 下坡减少大量参数
            lambda /= nu2;
        } else {
            // 上坡增加大量参数
            lambda *= nu1;
        }
        printf("\nlambda = %lf", lambda);
        printf("\nstep = %d", step);
        printf("\ne = %lf\n", e);
        step += 1;
    }

    free(J);
    free(r);
}


void trace(int n, double **A, double **trA) {
    // 构造对角线元素为A的对角线元素，其他元素均为0的矩阵
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j) {
                trA[i][j] = 0;
            } else {
                trA[i][j] = A[i][j];
            }
        }
    }
}


void output(int m, int n, double *beta) {
    int i, j;

    printf("beta = \n");
    for (i = 0; i < n; i++) {
        printf("%.10lf ", beta[i]);
    }
}