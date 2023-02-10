#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include "utils.h"

using namespace std;

double norm(double* vec, int sz){
	double ret = 0;
	for (int i = 0; i < sz; i++){
		ret += pow(vec[i], 2);
	}
	return sqrt(ret);
}

double* mvmulti(double** mat, double* vec, int sz){
	double* ret = (double*)malloc(sz * sizeof(double));
	for (int i = 0; i < sz; i ++){
		double sum = 0;
		for (int j = 0; j < sz; j++) sum += mat[i][j] * vec[j];
		ret[i] = sum;
	}
	return ret;
}

double* vvadd(double* v1, double v2, int sz, int sign){
	double* ret = (double*)malloc(sz * sizeof(double));
	for (int i = 0; i < sz; i ++) ret[i] = v1[i] + v2 * sign;
	return ret;
}

int main(int argc, char** argv){
	Timer t;
	long N = read_option<long>("-n", argc, argv);
	//long N = 100;
	long iter = 5000;
	double f = 1.0;
	double *u1 = (double *)malloc(N * sizeof(double));
	double *u2 = (double *)malloc(N * sizeof(double));

	for (int i = 0; i < N; i++) {
		u1[i] = 0;
		u2[i] = 0;
	}

	double** a = (double**)malloc(N * sizeof(double*));
	for (int i = 0; i < N; i ++){
		a[i] = (double *)malloc(N * sizeof(double));
		for (int j = 0; j < N; j++) a[i][j] = 0;
	}
	double d1 = -1 * pow(N+1,2);
	double d2 =  2 * pow(N+1, 2);
	for (int i = 1; i < N-1; i++){
		a[i][i-1] = d1;
		a[i][i]   = d2;
		a[i][i+1] = d1;
	}
	a[0][0] = d2;
	a[0][1] = d1;
	a[N-1][N-2] = d1;
	a[N-1][N-1] = d2;
	double* temp1 = mvmulti(a, u2, N);
	double* temp2 = vvadd(temp1, f, N, -1);
	free(temp1);
	double target = norm(temp2, N) / 1e4;
	free(temp2);
	int cur = 0;
	double cunorm = 0;
	double time = 0.0;

	cout << "================================================" << endl;
	cout << "Now using method 1: Jacobi's method." << endl;
	cout << "Iternation\t|Residual norm" << endl;
	t.tic();
	do{
		for (int i = 0; i < N; i++) u1[i] = u2[i];
		for (int i = 0; i < N; i++) {
			double tmp = 0.0;
			for (int j = 0; j < N; j++){
				if (j == i) continue;
				tmp+= a[i][j] * u1[j];
			}
			u2[i] = (f - tmp) / a[i][i];
		}
		cur++;
		double *temp1 = mvmulti(a, u2, N);
		double *temp2 = vvadd(temp1, f, N, -1);

		cunorm = norm(temp2, N);
		free(temp1);
		free(temp2);
		cout << cur << "\t\t|\t" << cunorm << endl;
	} while (cunorm > target && cur < iter);
	time = t.toc();
	cout << "Time spent =" << time << "s.\n" << endl;

	cunorm = 0;
	cout << "================================================" << endl;
	cout << "Now using method 2: Gauss-Seidel algo." << endl;
	cout << "Iternation\t|Residual norm" << endl;

	for (int i = 0; i < N; i++){
		u1[i] = 0;
		u2[i] = 0;
	}
	cur = 0;
	double *x;
	double *y;
	t.tic();
	do{
		x = (cur % 2 == 0) ? u1 : u2;
		y = (cur % 2 == 0) ? u2 : u1;
		for (int i = 0; i < N; i++){
			double tmp = f;
			for (int j = 0; j < i; j++) tmp -= a[i][j] * x[j];
			for (int j = i + 1; j < N; j++) tmp -= a[i][j] * y[j];
			tmp /= a[i][i];
			x[i] = tmp;
		}
		cur++;
		double *temp1 = mvmulti(a, x, N);
		double *temp2 = vvadd(temp1, f, N, -1);
		cunorm = norm(temp2, N);
		free(temp1);
		free(temp2);
		cout << cur << "\t\t|\t" << cunorm << endl;
	} while (cunorm > target && cur < iter);
	time = t.toc();
	cout << "Time spent =" << time << "s.\n" << endl;

	for (int i = 0; i < N; i++) free(a[i]);
	free(a);
	//free(f);
	free(u1);
	free(u2);
}
