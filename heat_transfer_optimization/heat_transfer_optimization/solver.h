#include <math.h>

void mult_mv(int n, double **mat, double *v, double *result){

    for (int i = 0; i < n; i++){
        result[i] = 0.;
        for (int j = 0; j < n; j++){
            result[i] += mat[i][j] * v[j];
        }
    }
}

double mult_vv(int n, double *v1, double *v2){

    double res = 0.;
    for (int i = 0; i < n ;i++){
        res += v1[i] * v2[i];
    }
    return res;
}

void mult_vs(int n, double *v, double s, double *res){

    for (int i = 0; i <n ; i++)
        res[i] = v[i] * s;
}

void sum_vv(int n, double *v1, double *v2, double *r){
    for (int i = 0; i < n; i++)
        r[i] = v1[i] + v2[i];
}

void sub_vv(int n, double *v1, double *v2, double *r){
    for (int i = 0; i < n; i++)
        r[i] = v1[i] - v2[i];
}

double norm_v(int n, double *v){

    return sqrt(mult_vv(n, v, v));
}

void residual(int n, double **A, double *x, double *b, double *r){

    mult_mv(n, A, x, r);
    sub_vv(n, b, r, r);
}

double residual_norm(int n, double **A, double *x, double *b){

    double result = 0., temp;
    for (int i = 0; i < n; i++)
    {
        temp = 0.;
        for (int j = 0; j < n; j++)
            temp += A[i][j] * x[j];

        temp -= b[i];
        result += temp * temp;
    }
    return sqrt(result);
}

void solver_CG(int n, double **A, double *b, double *x0, double *x){

    const int maxiter = 1000;
    const double maxres = 0.0001;
    double alpha = 0., beta = 0.;

    double *r = new double[n];
    double *r_old = new double[n];
    double *p = new double[n];
    double *Ap = new double[n];

    for (int i = 0; i < n; i++)
        x[i] = x0[i];

    residual(n, A, x, b, r);
    for (int j = 0; j < n; j++)
        p[j] = r[j];


    for (int i = 0; i < maxiter; i++)
    {
        const double residual_norm = norm_v(n, r);

        if (residual_norm < maxres)
        {
           // printf("CG solver converged in %d iterations with residual %.8f\n", i, residual_norm);
            delete[] r;
            delete[] r_old;
            delete[] p;
            delete[] Ap;
            return;
        }

        //printf("Current residual norm: %.8f\n", residual_norm);

        mult_mv(n, A, p, Ap);
        alpha = mult_vv(n, r, r)/mult_vv(n, p, Ap);

        for (int j = 0; j < n; j++){
            x[j] += alpha*p[j];
            r_old[j] = r[j];
            r[j] -= alpha*Ap[j];
        }

        beta = mult_vv(n, r, r) / mult_vv(n, r_old, r_old);

        for (int j = 0; j < n; j++)
            p[j] = r[j] + beta*p[j];

    }

    printf("CG solver failed to converge. The residual returned was %.8f\n", norm_v(n, r));
    delete[] r;
    delete[] r_old;
    delete[] p;
    delete[] Ap;
}