#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


double cal_error(int m, double *a, double *x, double *y) {
    // 计算误差
    int i, j;
    double e = 0.0, t;

    for (i = 0; i < m; i++) {
        t = a[0] + a[1] * sin(a[2] * x[i]) + a[3] * cos(a[4] * x[i]);
        e += pow(t - y[i], 2);
    }

    return sqrt(e);
}

void simulated_annealing(int n, int m, double *a, double *a_, double T, double *x, double *y) {
    /*
     模拟退火

     E 表示能量
     probability 表示退火概率
     choose 与probability表示保持状态概率
     */
    int i, j;
    double E, probability, choose;
    double alpha = 0.001;

    for (i = 0; i < n; i++) {
        a_[i] = ((rand() % ((int) (2 * alpha * 1e6 + 1))) + (a[i] - alpha) * 1e6) / 1e6;
    }
    if (cal_error(m, a_, x, y) <= cal_error(m, a, x, y)) {
        printf("new status:\n");
        for (i = 0; i < n; i++) {
            a[i] = a_[i];
            printf("%lf ", a[i]);
        }
        printf("\n");
    } else {
        E = cal_error(m, a, x, y) - cal_error(m, a_, x, y);
        probability = exp(E / T);
        choose = rand() % 1000 / 1000;
        if (choose <= probability) {
            printf("E<0, new status:\n");
            for (i = 0; i < n; i++) {
                a[i] = a_[i];
                printf("%lf ", a[i]);
            }
        }
        printf("\n");
    }
}


int main() {
    /*
     该程序用于使用模拟退火算法解决非线性最小二乘问题
     p(x)=a_0+a_1*sin(a_2*x)+a_3*cos(a_4*x)

     m 表示数据数目
     n 表示参数数目
     x 表示拟合点坐标
     y 表示拟合点函数值
     a 表示参数向量
     e 表示误差
     d 表示前后两次误差之差
     delta 表示d可允许最大值
     min_e 表示误差最小值
     T 表示初始温度
    */
    int m = 11, n = 5;
    double x[] = {-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    double y[] = {-0.8669, -0.2997, 0.3430, 1.0072, 1.6416, 2.2022, 2.6558, 2.9823, 3.1755, 3.2416, 3.1974};
    double a[] = {1.3, 2.2, 1.1, 0.9, 1.8};
    double *a_, *am;
    double eps = 0.001, e = 1, min_e = 1;
    double T = 5.0;
    int i, step = 0;
    srand((unsigned) time(NULL));

    a_ = (double *) malloc(n * sizeof(double));
    am = (double *) malloc(n * sizeof(double));

    while (e > eps) {
        step += 1;
        printf("step = %d\n", step);
        e = cal_error(m, a, x, y);
        if (e < min_e) {
            min_e = e;
            for (i = 0; i < n; i++) {
                am[i] = a[i];
            }
        }
        printf("error:\n%lf\n", e);
        simulated_annealing(n, m, a, a_, T, x, y);
        T *= 0.99;
        if (T < 1e-6) {
            printf("T = 0, can't find a answer");
            break;
        }
    }
    printf("\nmin_e: %lf\n", min_e);
    for (i = 0; i < n; i++) {
        printf("%lf ", am[i]);
    }

    free(a_);
    free(am);
}
