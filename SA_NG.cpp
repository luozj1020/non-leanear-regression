#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "vector_matrix.h"


void construct_J(int m, double *beta, double *x, double **J) {
    // 构造Jacobi矩阵J
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


void NGAlgorithm(int m, int n, double *beta, double alpha, double **J, double *r) {
    /*
     JT 为J的转置
     JTJ 为JT*J
     JTr 为JT*r
     */
    double **JT, **JTJ;
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
    JTr = (double *) malloc(n * sizeof(double));
    delta_beta = (double *) malloc(n * sizeof(double));

    matrix_tranform(m, n, J, JT);
    matrix_times_matrix(n, m, n, JT, J, JTJ);
    matrix_times_vector(n, m, JT, r, JTr);

    solve_linear(n, JTJ, JTr, delta_beta);
    num_times_vector(n, delta_beta, alpha, delta_beta);
    vector_add_vector(n, beta, delta_beta, beta);

    free(JT);
    free(JTJ);
    free(JTr);
    free(delta_beta);
}

double cal_error(int m, double *beta, double *x, double *y) {
    int i, j;
    double e = 0.0, t;

    for (i = 0; i < m; i++) {
        t = beta[0] + beta[1] * sin(beta[2] * x[i]) + beta[3] * cos(beta[4] * x[i]);
        e += pow(t - y[i], 2);
    }

    return sqrt(e);
}


void simulated_annealing(int n, int m, double *a, double T, double *x, double *y, double **J, double *r) {
    /*
     模拟退火

     E 表示能量
     probability 表示退火概率
     choose 与probability表示保持状态概率
     alpha 表示步长
     gamma 表示随机相邻状态的范围
     */
    int i, j;
    double E, probability, choose;
    double alpha = 0.1, gama = 0.5;
    double *_a;

    _a = (double *) malloc(n * sizeof(double));

    for (i = 0; i < n; i++) {
        _a[i] = a[i];
    }
    NGAlgorithm(m, n, a, alpha, J, r);
    if (cal_error(m, a, x, y) <= cal_error(m, _a, x, y)) {
        printf("new status:\n");
        for (i = 0; i < n; i++) {
            printf("%.10lf ", a[i]);
        }
        printf("\n");
    } else {
        E = cal_error(m, _a, x, y) - cal_error(m, a, x, y);
        probability = exp(E / T);
        choose = rand() % 1000 / 1000;
        if (choose > probability) {
            printf("E<0, new status:\n");
            for (i = 0; i < n; i++) {
                printf("%.10lf ", a[i]);
            }
        } else {
            printf("E<0, old status:\n");
            for (i = 0; i < n; i++) {
                a[i] = ((rand() % ((int) (2 * gama * 1e6 + 1))) + (_a[i] - gama) * 1e6) / 1e6;
                printf("%.10lf ", a[i]);
            }
        }
        printf("\n");
    }

    free(_a);
}


int main() {
    /*
     该程序用于解决非线性最小二乘问题
     p(x)=a_0+a_1*sin(a_2*x)+a_3*cos(a_4*x)

     m 表示数据数目
     n 表示参数数目
     x 表示拟合点坐标
     y 表示拟合点函数值
     beta 表示参数向量
     J 为Jacobi矩阵
     r 为残差
     alpha 表示梯度下降速率
     eps 表示可允许误差
     e 表示误差
     T 表示初始温度
    */
    int m = 11, n = 5;
    double x[] = {-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    double y[] = {-0.8669, -0.2997, 0.3430, 1.0072, 1.6416, 2.2022, 2.6558, 2.9823, 3.1755, 3.2416, 3.1974};
    // double beta[] = {0.1, 0.1, 0.1, 0.1, 0.1};
    double beta[] = {1.3, 2.2, 1.1, 0.9, 1.8};
    double **J;
    double *r;
    double eps = 1e-5, e = 10, T = 1.0;
    int step = 0;
    int i, j;
    srand((unsigned) time(NULL));

    // J是一个m*n矩阵
    J = (double **) malloc(sizeof(double *) * m);
    for (i = 0; i < m; i++) {
        J[i] = (double *) malloc(sizeof(double) * n);
    }
    r = (double *) malloc(m * sizeof(double));

    while (e > eps) {
        construct_J(m, beta, x, J);
        construct_r(m, x, y, beta, r);
        simulated_annealing(n, m, beta, T, x, y, J, r);
        e = cal_error(m, beta, x, y);
        printf("\ne = %.10lf\n", e);
        step += 1;
        printf("\nstep = %d\n", step);
        T *= 0.999;
        if (T < 1e-6) {
            break;
        }
    }

    free(J);
    free(r);
}
