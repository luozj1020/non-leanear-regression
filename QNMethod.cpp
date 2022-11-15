#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector_matrix.h"

double cal_error(int m, double *a, double *x, double *y);

double f(int i, double *a, double *x) {
    return a[0] + a[1] * sin(a[2] * x[i]) + a[3] * cos(a[4] * x[i]);
}

void nabla(int m, int n, double *a, double *x, double *y, double *nf) {
    // 计算f的梯度
    int i, j;

    for (i = 0; i < n; i++) {
        nf[i] = 0;
    }
    for (i = 0; i < m; i++) {
        nf[0] += 2 * (f(i, a, x) - y[i]);
        nf[1] += 2 * sin(a[2] * x[i]) * (f(i, a, x) - y[i]);
        nf[2] += 2 * a[1] * x[i] * cos(a[2] * x[i]) * (f(i, a, x) - y[i]);
        nf[3] += 2 * cos(a[4] * x[i]) * (f(i, a, x) - y[i]);
        nf[4] += (-2 * a[3] * x[i] * sin(a[4] * x[i]) * (f(i, a, x) - y[i]));
    }
}

double get_alpha(int m, int n, double *a, double *x, double *y, double *p) {
    // 寻找满足Wolfe条件的alpha
    double alpha = 1, c1 = 1e-4, c2 = 0.9;
    double l1 = 0, r1 = 0, l2 = 0, r2 = 0;
    double *ap, *xap, *nf, *nf_;
    int i, j;

    ap = (double *) malloc(sizeof(double) * n);
    xap = (double *) malloc(sizeof(double) * n);
    nf = (double *) malloc(sizeof(double) * n);
    nf_ = (double *) malloc(sizeof(double) * n);

    num_times_vector(n, p, alpha, ap);
    vector_add_vector(n, a, ap, xap);
    nabla(m, n, a, x, y, nf);
    nabla(m, n, xap, x, y, nf_);
    for (i = 0; i < m; i++) {
        l1 += pow(f(i, xap, x) - y[i], 2);
    }
    for (i = 0; i < m; i++) {
        r1 += pow(f(i, a, x) - y[i], 2);
    }
    r1 += c1 * vector_dot_vector(n, ap, nf);
    l2 = -vector_dot_vector(n, p, nf_);
    r2 = -c2 * vector_dot_vector(n, p, nf);

    while ((l1 - r1 > 1e-6) || (l2 - r2 > 1e-6)) {
        l1 = 0;
        r1 = 0;
        alpha *= 0.1;
        num_times_vector(n, p, alpha, ap);
        vector_add_vector(n, a, ap, xap);

        nabla(m, n, a, x, y, nf);
        nabla(m, n, xap, x, y, nf_);

        num_times_vector(n, p, alpha, ap);
        vector_add_vector(n, a, ap, xap);
        nabla(m, n, a, x, y, nf);
        nabla(m, n, xap, x, y, nf_);
        for (i = 0; i < m; i++) {
            l1 += pow(f(i, xap, x) - y[i], 2);
        }
        for (i = 0; i < m; i++) {
            r1 += pow(f(i, a, x) - y[i], 2);
        }
        r1 += c1 * vector_dot_vector(n, ap, nf);
        l2 = -vector_dot_vector(n, p, nf_);
        r2 = -c2 * vector_dot_vector(n, p, nf);
    }

    free(ap);
    free(xap);
    free(nf);
    free(nf_);

    return alpha;
}

void QNMethod(int m, int n, double beta, double *a, double *x, double *y, double eps) {
    /*
     nf 表示f(x_k)
     nf_ 表示f(x_{k+1})
     p 表示方向向量
     s 表示a_k*p_k
     a_ 表示x_{k+1}
     yk 表示y_k=\nabla f(x_{k+1})-\nabla f(x_{k})
     ss 表示s_k*s_k^T
     Bys 表示B_k*y_k
     syB 表示s_k*y_k^T*B_k
     yB 表示y_k^T*B_k
     By 表示B_k*y_k
     B 表示B_k
     B_ 表示B_{k+1}
     alpha 表示\alpha_k
     e 表示误差
     */
    double *nf, *nf_, *p, *s, *a_, *yk;
    double **ss, **Bys, **syB;
    double *yB, *By;
    double sy, yBy;
    double **B, **B_;
    double alpha, e = 1.0;
    int i, j, step = 0;

    B = (double **) malloc(sizeof(double *) * n);
    for (i = 0; i < n; i++) {
        B[i] = (double *) malloc(sizeof(double) * n);
    }
    B_ = (double **) malloc(sizeof(double *) * n);
    for (i = 0; i < n; i++) {
        B_[i] = (double *) malloc(sizeof(double) * n);
    }
    ss = (double **) malloc(sizeof(double *) * n);
    for (i = 0; i < n; i++) {
        ss[i] = (double *) malloc(sizeof(double) * n);
    }
    Bys = (double **) malloc(sizeof(double *) * n);
    for (i = 0; i < n; i++) {
        Bys[i] = (double *) malloc(sizeof(double) * n);
    }
    syB = (double **) malloc(sizeof(double *) * n);
    for (i = 0; i < n; i++) {
        syB[i] = (double *) malloc(sizeof(double) * n);
    }
    nf = (double *) malloc(sizeof(double) * n);
    nf_ = (double *) malloc(sizeof(double) * n);
    p = (double *) malloc(sizeof(double) * n);
    s = (double *) malloc(sizeof(double) * n);
    a_ = (double *) malloc(sizeof(double) * n);
    yk = (double *) malloc(sizeof(double) * n);
    yB = (double *) malloc(sizeof(double) * n);
    By = (double *) malloc(sizeof(double) * n);

    for (i = 0; i < n; i++) {
        B[i][i] = beta;
    }

    while (e > eps) {
        step += 1;
        printf("step = %d\n", step);
        // step1: 计算p_k
        nabla(m, n, a, x, y, nf);
        matrix_times_vector(n, n, B, nf, p);
        num_times_vector(n, p, -1.0, p);
        // step2: 计算alpha
        alpha = get_alpha(m, n, a, x, y, p);
        printf("alpha = %lf\n", alpha);
        // step3: 计算x_{k+1}
        num_times_vector(n, p, alpha, s);
        vector_add_vector(n, a, s, a_);
        printf("a_ = \n");
        for (i = 0; i < n; i++) {
            printf("%lf ", a_[i]);
        }
        printf("\n");
        e = cal_error(m, a_, x, y);
        printf("e = %lf\n", e);
        // step4: 计算y_k
        nabla(m, n, a_, x, y, nf_);
        vector_minus_vector(n, nf_, nf, yk);
        // step5: 计算B_{k+1}
        sy = vector_dot_vector(n, s, yk);
        matrix_times_vector(n, n, B, yk, By);
        yBy = vector_dot_vector(n, yk, By);
        vector_times_vector(n, s, s, ss);
        vector_times_vector(n, By, s, Bys);
        vector_times_matrix(n, n, B, yk, yB);
        vector_times_vector(n, s, yB, syB);
        // 计算第二项
        num_times_matrix(n, n, ss, sy + yBy, ss);
        matrix_divide_num(n, n, ss, pow(sy, 2), ss);
        // 计算第三项
        matrix_add_matrix(n, n, Bys, syB, syB);
        matrix_divide_num(n, n, syB, sy, syB);
        // 计算B_{k+1}
        matrix_add_matrix(n, n, B, ss, B_);
        matrix_minus_matrix(n, n, B_, syB, B_);
        // 更新数据
        copy_vector(n, a_, a);
        copy_matrix(n, n, B_, B);
    }


    free(B);
    free(B_);
    free(ss);
    free(Bys);
    free(syB);
    free(nf);
    free(nf_);
    free(p);
    free(s);
    free(a_);
    free(yk);
    free(yB);
    free(By);
}


int main() {
    /*
     该程序用于使用准牛顿法解决非线性最小二乘问题
     p(x)=a_0+a_1*sin(a_2*x)+a_3*cos(a_4*x)

     m 表示数据数目
     n 表示参数数目
     x 表示拟合点坐标
     y 表示拟合点函数值
     a 表示参数向量
     eps 表示可允许误差
     */
    int m = 11, n = 5;
    double x[] = {-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    double y[] = {-0.8669, -0.2997, 0.3430, 1.0072, 1.6416, 2.2022, 2.6558, 2.9823, 3.1755, 3.2416, 3.1974};
    double a[] = {0.1, 0.1, 0.1, 0.1, 0.1};
    // double a[] = {1.3, 2.2, 1.1, 0.9, 1.8};
    double beta = 1.0, eps = 0.00001;

    QNMethod(m, n, beta, a, x, y, eps);
}


double cal_error(int m, double *a, double *x, double *y) {
    int i, j;
    double e = 0.0, t;

    for (i = 0; i < m; i++) {
        t = a[0] + a[1] * sin(a[2] * x[i]) + a[3] * cos(a[4] * x[i]);
        e += pow(t - y[i], 2);
    }

    return sqrt(e);
}