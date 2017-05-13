#include <stdio.h>

#define MAX_SIZE 10
#include <math.h>
void swap(double *x, double *y) {
	double t = *x;
	*x = *y;
	*y = t;
}
int Solve(int n, double a[][MAX_SIZE], double x[], double b[]) {
	double y[MAX_SIZE];
	for(int i = 0; i < n; i++) {
		double sum = 0;
		for(int j = 0; j < i; j++)
			sum += a[i][j] * y[j];
		y[i] = b[i] - sum;
	}
	for(int i = n-1; i >= 0; i--) {
		double sum = 0;
		for(int j = n-1; j > i; j--)
			sum += a[i][j] * x[j];
		x[i] = (y[i] - sum) / a[i][i];
	}
}
int doolittle_LU(int n, double a[][MAX_SIZE]) {
	for(int i = 0; i < n; i++) {
		if(fabs(a[i][i]) < 1e-15)
			return 0;
		for(int j = i; j < n; j++) {
			double sum = 0;
			for(int k = 0; k < i; k++)
				sum += a[i][k] * a[k][j];
			a[i][j] -= sum;
		}
		for(int j = i+1; j < n; j++) {
			double sum = 0;
			for(int k = 0; k < i; k++) 
				sum += a[j][k] * a[k][i];
			a[j][i] = (a[j][i] - sum) / a[i][i];
		}
	}
	return 1;
}
int EigenV(int n, double b[][MAX_SIZE], double *lambda, double v[], double TOL, int MAXN) {
	double l0 = *lambda;
	double a[MAX_SIZE][MAX_SIZE];
	int k = 0;
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++) a[i][j] = b[i][j];
	for(int i = 0; i < n; i++)
		a[i][i] -= *lambda;
	if(!doolittle_LU(n, a)) {
		return -1;
	}
	for (int i = 0; i < n; i++) {
		if (a[i][i] == 0) {
			return -1;
		}
	}
	for(int i = 0; i < n; i++)
		if(fabs(v[i]) > fabs(v[k])) k = i;
	double tmp = v[k];
	for(int i = 0; i < n; i++)
		v[i] /= tmp;
	double err = 0;
	double u[MAX_SIZE];
	for(int times = 0; times < MAXN; times++) {
		for(int i = 0; i < n; i++)
			u[i] = v[i];
		// for(int i = 0; i < n; i++)
		// 	printf("%lf ", u[i]);
		// puts("----");
		// for(int i = 0; i < n; i++)
		// 	printf("%lf ", v[i]);
		// puts("===");
		Solve(n, a, v, u);
		*lambda = 1.0 / v[k];

		for(int i = 0; i < n; i++)
			if(fabs(v[i]) > fabs(v[k])) k = i;
		double tmp = v[k];
		for(int i = 0; i < n; i++)
			v[i] /= tmp;
		
		err = 0;
		for(int i = 0; i < n; i++) {
			if(fabs(u[i] - v[i]) > err) err = fabs(u[i] - v[i]);
		}
		// printf("%.12lf\n", err);
		// printf("%.12lf\n", TOL);
		if (err < TOL) {
			*lambda += l0;
			return 1;
		}
		// printf("%d\n", times);
	}
	return 0;
}
int main()
{
	// test();
	// return 0;
	int n, MAXN, m, i, j, k;
	double a[MAX_SIZE][MAX_SIZE], v[MAX_SIZE];
	double lambda, TOL;

	scanf("%d", &n);
	for (i=0; i<n; i++) 
		for (j=0; j<n; j++) 
			scanf("%lf", &a[i][j]);
	scanf("%lf %d", &TOL, &MAXN);
	scanf("%d", &m);
	for (i=0; i<m; i++) {
		scanf("%lf", &lambda);
		for (j=0; j<n; j++)
			scanf("%lf", &v[j]);
		switch (EigenV(n, a, &lambda, v, TOL, MAXN)) {
			case -1: 
				printf("%12.8f is an eigenvalue.\n", lambda );
				break;
			case 0:
				printf("Maximum number of iterations exceeded.\n");
				break;
			case 1:
				printf("%12.8f\n", lambda );
				for (k=0; k<n; k++)
					printf("%12.8f ", v[k]);
				printf("\n");
				break;
		}
	}

	return 0;
}

/* Your function will be put here */