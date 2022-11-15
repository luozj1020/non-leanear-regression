#include <stdio.h>
#include <math.h>

void matrix_tranform(int m, int n, double **A, double **B) {
    // 矩阵转置
    int i, j;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            B[j][i] = A[i][j];
        }
    }
}

void copy_matrix(int m, int n, double **A, double **B) {
    // 矩阵复制
    int i, j;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            B[i][j] = A[i][j];
        }
    }
}

void matrix_add_matrix(int m, int n, double **A, double **B, double **C) {
    // 矩阵相加
    int i, j;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

void matrix_minus_matrix(int m, int n, double **A, double **B, double **C) {
    // 矩阵相减
    int i, j;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}

void matrix_times_matrix(int m, int n, int p, double **A, double **B, double **C) {
    // 矩阵相乘
    int i, j, k;

    for (i = 0; i < m; i++) {
        for (j = 0; j < p; j++) {
            C[i][j] = 0.0;
            for (k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void num_times_matrix(int m, int n, double **A, double b, double **C) {
    // 数字乘矩阵
    int i, j;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            C[i][j] = b * A[i][j];
        }
    }
}

void matrix_divide_num(int m, int n, double **A, double b, double **C) {
    // 矩阵除数字
    int i, j;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            C[i][j] = A[i][j] / b;
        }
    }
}

void matrix_times_vector(int m, int n, double **A, double *b, double *c) {
    // 矩阵乘向量
    int i, j;

    for (i = 0; i < m; i++) {
        c[i] = 0.0;
        for (j = 0; j < n; j++) {
            c[i] += A[i][j] * b[j];
        }
    }
}

void vector_times_matrix(int m, int n, double **A, double *b, double *c) {
    // 向量乘矩阵
    int i, j;

    for (i = 0; i < n; i++) {
        c[i] = 0;
        for (j = 0; j < m; j++) {
            c[i] += b[j] * A[j][i];
        }
    }
}

void copy_vector(int n, double *a, double *b) {
    // 向量复制
    int i;

    for (i = 0; i < n; i++) {
        b[i] = a[i];
    }
}

void vector_add_vector(int n, double *a, double *b, double *c) {
    // 向量相加
    int i;

    for (i = 0; i < n; i++) {
        c[i] = a[i] + b[i];
    }
}

void vector_minus_vector(int n, double *a, double *b, double *c) {
    // 向量相减
    int i;

    for (i = 0; i < n; i++) {
        c[i] = a[i] - b[i];
    }
}

void num_times_vector(int n, double *a, double b, double *c) {
    // 数字乘向量
    int i;

    for (i = 0; i < n; i++) {
        c[i] = b * a[i];
    }
}

void vector_divide_num(int n, double *a, double b, double *c) {
    // 向量除数字
    int i;

    for (i = 0; i < n; i++) {
        c[i] = a[i] / b;
    }
}

double vector_dot_vector(int n, double *a, double *b) {
    // 向量内积，即(1*n)*(n*1)
    int i;
    double t = 0;

    for (i = 0; i < n; i++) {
        t += a[i] * b[i];
    }

    return t;
}

void vector_times_vector(int n, double *a, double *b, double **M) {
    // 向量相乘，但是(n*1)*(1*n)
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            M[i][j] = a[i] * b[j];
        }
    }
}

void solve_linear(int n, double **A, double *b, double *x) {
    // 求解线性方程组Ax=b
    int i, j, k, r, kk, ik, ki, ri;
    double iMax, t, Aii;

    for (k = 0; k < n - 1; k++) {
        r = k;
        kk = k * n + k;
        iMax = fabs(A[k][k]);
        for (i = k + 1; i < n; i++) {
            t = fabs(A[i][k]);
            if (t > iMax) {
                r = i;
                iMax = t;
            }
        }   // 选列主元
        if (r != k) {
            for (i = k; i < n; i++) {
                ki = k * n + i;
                ri = r * n + i;
                t = A[k][i];
                A[k][i] = A[r][i];
                A[r][i] = t;
            } // 交换矩阵 A 的 r,k 两行元素
            t = b[k];
            b[k] = b[r];
            b[r] = t;  // 交换 b 的 r,k 两行元素
        }

        if (fabs(A[k][k]) < 1e-12) {
            printf("fail\n");
            return;
        }
        for (i = k + 1; i < n; i++) {
            ik = i * n + k;
            A[i][k] /= A[k][k];
            b[i] -= A[i][k] * b[k];
            for (j = k + 1; j < n; j++) {
                A[i][j] -= A[i][k] * A[k][j];
            }
            A[i][k] = 0.0;
        }
    }

    kk = k * n + k;
    if (fabs(A[k][k]) < 1e-12) {
        printf("fail\n");
        return;
    }

    Aii = A[n - 1][n - 1];
    if (fabs(Aii) < 1e-12) {
        printf("fail\n");
        return;
    } else {
        x[n - 1] = b[n - 1];
        if (Aii != 1.0) {
            x[n - 1] /= Aii;
        }
    }

    for (i = n - 2; i >= 0; i--) {
        Aii = A[i][i];
        if (fabs(Aii) < 1e-12) {
            printf("fail\n");
            return;
        } else {
            x[i] = 0.0;
            for (j = i + 1; j < n; j++) {
                x[i] += A[i][j] * x[j];
            }
            x[i] = b[i] - x[i];
            if (Aii != 1.0) {
                x[i] /= Aii;
            }
        }
    }
}