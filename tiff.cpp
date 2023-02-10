#include <stdio.h>
#include "utils.h"
#include <iostream>
#include <numeric>
#include <cmath>

using namespace std;

void jacobi(long N, double ite, double f)
{
    double h = 1.0 / (N + 1);                            // length of step
    double *u_kk = (double *)malloc(N * sizeof(double)); // store the vector u_{k+1}
    double *u_k = (double *)malloc(N * sizeof(double));  // store the vector u_{k}
    double *Au = (double *)malloc(N * sizeof(double));
    double err;
    double A_off = -1 / pow(h, 2);
    double A_diag = 2 / pow(h, 2);
    double sum = 0;
    double err0 = sqrt(N);
    memset(u_k, 0, sizeof(*u_k));
    memset(u_kk, 0, sizeof(*u_k));

    for (int k = 0; k < ite; k++)
    {
        err = 0;

        for (int i = 0; i < N; i++)
        {

            if (i == 0)
            {
                sum = A_off * u_k[i + 1];
            }
            if (i == N - 1)
            {
                sum = A_off * u_k[i - 1];
            }
            if (i != 0 && i != N - 1)
                sum = (A_off * u_k[i - 1] + A_off * u_k[i + 1]);
            u_kk[i] = 1 / A_diag * (f - sum); // jocobi iteration
        }

        for (int i = 0; i < N; i++)
        {
            u_k[i] = u_kk[i];
            if (i > 0 && i < N - 1)
            {
                Au[i] = A_off * u_kk[i - 1] + A_diag * u_kk[i] + A_off * u_kk[i + 1];
            }
            if (i == 0)
            {
                Au[i] = A_diag * u_kk[i] + A_off * u_kk[i + 1];
            }
            if (i == N - 1)
            {
                Au[i] = A_off * u_kk[i - 1] + A_diag * u_kk[i];
            }
            err = err + pow(Au[i] - 1, 2);
        }

        err = sqrt(err);

        if (err < err0 / pow(10, 4))
        {
            cout << "number of iteration=" << k + 1 << endl
                 << "residual=" << err << endl;
            break;
        }
        if (k == ite - 1)
        {
            cout << "iteation finished for k=5000" << endl
                 << "residual=" << err << endl;
        }
    }

    free(u_k);
    free(u_kk);
    free(Au);
}

void gauss_seidal(long N, double ite, double f)
{
    double h = 1.0 / (N + 1);                           // length of step
    double *u_k = (double *)malloc(N * sizeof(double)); // store the vector u_{k}
    double *Au = (double *)malloc(N * sizeof(double));
    double err;
    double A_off = -1 / pow(h, 2);
    double A_diag = 2 / pow(h, 2);
    double sum = 0;
    double err0 = sqrt(N);

    memset(u_k, 0, sizeof(*u_k));

    for (int k = 0; k < ite; k++)
    {
        err = 0;

        for (int i = 0; i < N; i++)
        {
            if (i == 0)
            {
                sum = A_off * u_k[i + 1];
            }
            if (i == N - 1)
            {
                sum = A_off * u_k[i - 1];
            }
            if (i != 0 && i != N - 1)
                sum = (A_off * u_k[i - 1] + A_off * u_k[i + 1]);
            u_k[i] = 1 / A_diag * (f - sum); // gauss seidal iteration
        }

        for (int i = 0; i < N; i++)
        {
            if (i > 0 && i < N - 1)
            {
                Au[i] = A_off * u_k[i - 1] + A_diag * u_k[i] + A_off * u_k[i + 1];
            }
            if (i == 0)
            {
                Au[i] = A_diag * u_k[i] + A_off * u_k[i + 1];
            }
            if (i == N - 1)
            {
                Au[i] = A_off * u_k[i - 1] + A_diag * u_k[i];
            }
            err = err + pow(Au[i] - 1, 2);
        }

        err = sqrt(err);

        if (err < err0 / pow(10, 4))
        {
            cout << "iteration time=" << k + 1 << endl
                 << "residual=" << err << endl;
            break;
        }
        if (k == ite - 1)
        {
            cout << "iteation finished for k=" << ite << endl
                 << "residual=" << err << endl;
        }
    }

    free(u_k);
    free(Au);
}

int main()
{

    long N = 100000; // demension of matrix A
    int ite = 5000;  // number of iteration
    double f = 1;
    double time;

    cout << "when N =" << N << endl;
    Timer t;
    t.tic();
    cout << "-----------" << endl
         << "Jacobi Method:" << endl;
    jacobi(N, ite, f);
    time = t.toc();
    cout << "time = " << time << endl;

    t.tic();
    cout << "-----------" << endl
         << "Gauss-Seidal Method:" << endl;
    gauss_seidal(N, ite, f);
    time = t.toc();
    cout << "time = " << time << endl;

    return 0;
}