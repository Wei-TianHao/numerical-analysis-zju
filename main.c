#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#pragma warning(disable : 4996)
#define MAX_SIZE 10

int EigenV(int n, double a[][MAX_SIZE], double *lambda, double v[], double TOL, int MAXN);

int main()
{
	int n, MAXN, m, i, j, k;
	double a[MAX_SIZE][MAX_SIZE], v[MAX_SIZE];
	double lambda, TOL;

	scanf("%d", &n);
	for (i = 0; i<n; i++)
		for (j = 0; j<n; j++)
			scanf("%lf", &a[i][j]);
	scanf("%lf %d", &TOL, &MAXN);
	scanf("%d", &m);
	for (i = 0; i<m; i++) {
		scanf("%lf", &lambda);
		for (j = 0; j<n; j++)
			scanf("%lf", &v[j]);
		switch (EigenV(n, a, &lambda, v, TOL, MAXN)) {
		case -1:
			printf("%12.8f is an eigenvalue.\n", lambda);
			break;
		case 0:
			printf("Maximum number of iterations exceeded.\n");
			break;
		case 1:
			printf("%12.8f\n", lambda);
			for (k = 0; k<n; k++)
				printf("%12.8f ", v[k]);
			printf("\n");
			break;
		}
	}
	system("pause");
	return 0;
}
void getCofactor(double A[MAX_SIZE][MAX_SIZE], double temp[MAX_SIZE][MAX_SIZE], int p, int q, int n) {
	int i = 0;
	int j = 0;
	for (int row = 0; row < n; row++) {
		for (int col = 0; col < n; col++) {
			if (row != p && col != q) {
				temp[i][j++] = A[row][col];
				if (j == n - 1) {
					j = 0;
					i++;
				}
			}
		}
	}
}
double getDeterminent(double a[][MAX_SIZE], int n) {
	if (n == 1) {
		return a[0][0];
	}
	if (n == 2) {
		return a[0][0] * a[1][1] - a[0][1] * a[1][0];
	}
	double ret = 0.;
	int sign = 1;
	double temp[MAX_SIZE][MAX_SIZE];
	for (int i = 0; i < n; i++) {
		getCofactor(a, temp, 0, i, n);
		ret += sign * a[0][i] * getDeterminent(temp, n - 1);
		sign = -sign;
	}
	return ret;
}
void getAStar(double a[][MAX_SIZE], double ret[][MAX_SIZE], int n) {
	if (n == 1) {
		ret[0][0] = 1.;
	}
	int sign = 1;
	double temp[MAX_SIZE][MAX_SIZE];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			getCofactor(a, temp, i, j, n);
			sign = ((i + j) % 2 == 0) ? 1 : -1;
			ret[j][i] = sign * (getDeterminent(temp, n - 1));
		}
	}
}
int inverse(double a[MAX_SIZE][MAX_SIZE], double ret[MAX_SIZE][MAX_SIZE], int n) {
	double det = getDeterminent(a, n);
	if (det == 0) {
		return 0;
	}
	double adj[MAX_SIZE][MAX_SIZE];
	getAStar(a, adj, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			ret[i][j] = adj[i][j] / det;
		}
	}
	return 1;
}
void multiply(double a[][MAX_SIZE],double x[MAX_SIZE], double result[MAX_SIZE], int n) {
	for (int i = 0; i < n; i++) {
		result[i] = 0;
		for (int j = 0; j < n; j++) {
			result[i] += a[i][j] * x[j];
		}
	}
}
double normalize(double a[MAX_SIZE], int n) {
	double max = 0;
	int sign = 1;
	for (int i = 0; i < n; i++) {
		if (fabs(a[i]) > max) {
			max = fabs(a[i]);
			if (a[i] < 0) {
				sign = -1;
			}
			else {
				sign = 1;
			}
		}
	}
	max = max * sign;
	for (int i = 0; i < n; i++) {
		a[i] /= max;
	}
	return max;
}
int dooLittle(int n, double a[MAX_SIZE][MAX_SIZE]) {
	double sum = 0;
	for (int i = 0; i < n; i++) {
		if (a[i][i] == 0) {
			return 0;
		}
		for (int j = i; j < n; j++) {
			sum = 0;
			for (int k = 0; k < i; k++) {
				sum += a[i][k] * a[k][j];
			}
			a[i][j] -= sum;
		}
		for (int j = i + 1; j < n; j++) {
			sum = 0;
			for (int k = 0; k < i; k++) {
				sum += a[j][k] * a[k][i];
			}
			a[j][i] = (a[j][i] - sum) / a[i][i];
		}
		
	}
	return 1;
}
int Solve(int n, double a[][MAX_SIZE], double v[], double b[]) {
	for (int i = 0; i < n; i++) {
		double sum = 0;
		for (int j = 0; j < i; j++)
			sum += a[i][j] * v[j];
		v[i] -= sum;
	}
	for (int i = n - 1; i >= 0; i--) {
		double sum = 0.0;
		for (int j = i + 1; j < n; j++)
			sum += a[i][j] * v[j];
		v[i] = (v[i] - sum) / a[i][i];
	}
}
int EigenV(int n, double b[][MAX_SIZE], double *lambda, double v[], double TOL, int MAXN) {
	double a[MAX_SIZE][MAX_SIZE];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			a[i][j] = b[i][j];
		}
	}
	for (int i = 0; i < n; i++) {
		a[i][i] -= *lambda;
	}
	if (!dooLittle(n, a))
		return -1;
	for (int i = 0; i < n; i++) {
		if (a[i][i] == 0) {
			return -1;
		}
	}
	double max = 0;
	int k = 0;
	for (int i = 0; i < n; i++) {
		if (fabs(v[i]) > fabs(max)) {
			max = v[i];
			k = i;
		}
	}
	for (int i = 0; i < n; i++) {
		v[i] /= max;
	}
	double u[MAX_SIZE];
	for (int k = 1; k <= MAXN; k++) {
		for (int i = 0; i < n; i++) {
			u[i] = v[i];
		}

		Solve(n, a, v, u);

		double mu = v[k];
		max = 0;
		for (int i = 0; i < n; i++) {
			if (fabs(v[i]) > fabs(max)) {
				max = v[i];
				k = i;
			}
		}
		for (int i = 0; i < n; i++) {
			v[i] /= max;
		}
		double error = 0;
		for (int i = 0; i < n; i++) {
			if (fabs(u[i] - v[i]) > error) {
				error = fabs(u[i] - v[i]);
			}
		}
		if (error < TOL) {
			*lambda = 1. / mu + *lambda;
			return 1;
		}
	}
	return 0;
}