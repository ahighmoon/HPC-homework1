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
		for (int j = 0; j < sz; j++){
			sum += mat[i][j] * vec[j];
		}
		ret[i] = sum;
	}
	return ret;
}

double* vvadd(double* v1, double* v2, int sz){
	double* ret = (double*)malloc(sz * sizeof(double));
	for (int i = 0; i < sz; i ++){
		ret[i] = v1[i] + v2[i];
	}
	return ret;
}

int main(int argc, char** argv){
	//Timer t;
	long N = read_option<long>("-n", argc, argv);
	double h = 1.0/ (N+1);
	long iter = 5000;
	double* f = (double*)malloc(N * sizeof(double));
	double* u = (double*)malloc(N * sizeof(double));
	for (int i = 0; i < N; i++) {
		f[i] = 1;
		u[i] = 0;
	}
	double** a = (double**)malloc(N * N * sizeof(double));
	memset(a, 0, N * N);
	double d1 = -(N+1)^2;
	double d2 = 2 * (N+1)^2;
	for (int i = 1; i < N-1; i++){
		a[i][i-1] = d1;
		a[i][i]   = d2;
		a[i][i+1] = d1;
	}
	a[0][0] = d2;
	a[0][1] = d1;
	a[N-1][N-2] = d1;
	a[N-1][N-1] = d2;
	double* temp1 = mvmulti(a, u, N);
	double* temp2 = vvadd(temp1, f, N);
	free(temp1);
	double target = norm(temp2, N) / 1e4;
	free(temp2);
	//for (int it = 0; it < iter; it++){
	int cur = 0;
	double curnorm = 0;
	do{
		for (int i = 0; i < N; i++){
			double tmp = 0.0;
			for (int j = 0; j < N; j++){
				if (j == i) continue;
				tmp+= a[i][j] * u[j];
			}
			u[i] = (f[i] - tmp) / a[i][i];
		}
		cur++;
		curnorm = norm(u, N);
		cout << "iter = " << cur+1 << ", norm of the residual = " << curnorm << endl;
	} while (curnorm > target && cur < iter);

	free(a);
	free(f);
	free(u);

}
