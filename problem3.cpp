#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include "utils.h"

using namespace std;

double norm(double* vec, int sz){
	double ret = 0;
	for (int i = 0; i < sz; i++) ret += pow(vec[i], 2);
	return sqrt(ret);
}

double* mvmulti(double* vec, double d1, double d2, int sz){
	double* ret = (double*)malloc(sz * sizeof(double));
	for (int i = 0; i < sz; i ++){
		double sum = d2 * vec[i];
		if (i == 0) sum += d1 * vec[i+1];
		else if (i == sz - 1) sum += d1 * vec[i-1];
		else sum += d1 * vec[i - 1] + d1 * vec[i + 1];
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
	//long N = read_option<long>("-n", argc, argv);
	long N = 100000;
	double h = 1.0 / (N + 1);
	long iter = 5000;
	double f = 1.0;
	double *u1 = (double *)malloc(N * sizeof(double));
	double *u2 = (double *)malloc(N * sizeof(double));

	for (int i = 0; i < N; i++) {
		u1[i] = 0;
		u2[i] = 0;
	}
	double d1 = -1 / pow(h, 2);
	double d2 = 2 / pow(h, 2);
	//double d1 = -1 * pow(N+1, 2);
	//double d2 =  2 * pow(N+1, 2);
	double* temp1 = mvmulti(u2, d1, d2, N);
	double* temp2 = vvadd(temp1, f, N, -1);
	free(temp1);
	double target = norm(temp2, N) / 1e4;
	free(temp2);
	int cur = 0;
	double cunorm = 0;
	double time = 0.0;

	//cout << "================================================" << endl;
	cout << "Now using method 1: Jacobi's method." << endl;
	//cout << "Iternation\t|Residual norm" << endl;
	t.tic();
	do{
		for (int i = 0; i < N; i++) u1[i] = u2[i];
		for (int i = 0; i < N; i++) {
			double tmp = 0;
			if (i == 0) {tmp = d1 * u1[1];}
			else if (i == N-1) {tmp = d1 * u1[N-2];}
			else {tmp = d1 * u1[i - 1] + d1 * u1[i + 1];}
			u2[i] = 1 / d2 * (f - tmp);
		}
		cur++;
		double *temp1 = mvmulti(u2, d1, d2, N);
		double *temp2 = vvadd(temp1, f, N, -1);
		cunorm = norm(temp2, N);
		free(temp1);
		free(temp2);
		cout << cur << "\t\t|\t" << cunorm << endl;
	} while (cunorm > target && cur < iter);
	cout << cur << "\t\t|\t" << cunorm << endl;
	time = t.toc();
	cout << "Time spent =" << time << "s.\n" << endl;

	cunorm = 0;
	//cout << "================================================" << endl;
	cout << "Now using method 2: Gauss-Seidel algo." << endl;
	//cout << "Iternation\t|Residual norm" << endl;

	for (int i = 0; i < N; i++) u1[i] = 0;
	cur = 0;
	t.tic();
	do{
		for (int i = 0; i < N; i++){
			double tmp = f;
			if (i == 0)
				tmp -= d1 * u1[1];
			else if (i == N - 1)
				tmp -= d1 * u1[N - 2];
			else
				tmp -= (d1 * u1[i - 1] + d1 * u1[i + 1]);
			u1[i] = tmp / d2;
		}
		cur++;
		double *temp1 = mvmulti(u1, d1, d2, N);
		double *temp2 = vvadd(temp1, f, N, -1);
		cunorm = norm(temp2, N);
		free(temp1);
		free(temp2);
		//cout << cur << "\t\t|\t" << cunorm << endl;
	} while (cunorm > target && cur < iter);
	time = t.toc();
	cout << cur << "\t\t|\t" << cunorm << endl;
	cout << "Time spent =" << time << "s.\n" << endl;

	free(u1);
	free(u2);
}
