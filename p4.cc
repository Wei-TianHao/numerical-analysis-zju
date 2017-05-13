#include <stdio.h>
#include <math.h>

#define MAX_SIZE 10
#define bound pow(2, 127)
#define ZERO 1e-9 /* X is considered to be 0 if |X|<ZERO */

int Jacobi( int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN ) {
	double diff, bd;
	double nx[MAX_SIZE];
	for(int i = 0; i < n; i++) {
		int k = i;
		for(int j = i; j < n; j++) {
			if(fabs(a[j][i]) > fabs(a[k][i])) k = j;
		}
		if(fabs(a[k][i]) < ZERO) {
			k = i;
			for(int j = i; j >=0; j--) {
				if(fabs(a[j][i]) > fabs(a[k][i])) k = j;
			}
			if(fabs(a[k][i]) < ZERO) {
				return -1;
			}
			for(int j = 0; j < n; j++) a[i][j] += a[k][j];
			b[i] += b[k];
		}
		else {
			for(int j = 0; j < n; j++) {
				int tmp = a[i][j];
				a[i][j] = a[k][j];
				a[k][j] = tmp;
			}
			int tmp = b[i];
			b[i] = b[k];
			b[k] = tmp;
		}
	}
	for(int it = 0; it < MAXN; it++){
		for(int i = 0; i < n; i++) {
			nx[i] = 0;
			for(int j = 0; j < n; j++) if(i != j) {
				nx[i] -= a[i][j] * x[j];
			}
			nx[i] += b[i];
			nx[i] /= a[i][i];
		}
		diff = bd = 0;
		for(int i = 0; i < n; i++){
			if(fabs(nx[i] - x[i]) > diff) diff = fabs(nx[i] - x[i]);
			if(fabs(nx[i]) > bd) bd = fabs(nx[i]);
			x[i] = nx[i];
		}
		if(bd > bound) return -2;
		if(diff < TOL) return it+1;
	}
	return 0;
}

int Gauss_Seidel( int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN ) {
	double diff, bd;
	for(int i = 0; i < n; i++) {
		int k = i;
		for(int j = i; j < n; j++) {
			if(fabs(a[j][i]) > fabs(a[k][i])) k = j;
		}
		if(fabs(a[k][i]) < ZERO) {
			k = i;
			for(int j = i; j >=0; j--) {
				if(fabs(a[j][i]) > fabs(a[k][i])) k = j;
			}
			if(fabs(a[k][i]) < ZERO) {
				return -1;
			}
			for(int j = 0; j < n; j++) a[i][j] += a[k][j];
			b[i] += b[k];
		}
		else {
			for(int j = 0; j < n; j++) {
				int tmp = a[i][j];
				a[i][j] = a[k][j];
				a[k][j] = tmp;
			}
			int tmp = b[i];
			b[i] = b[k];
			b[k] = tmp;
		}
	}
	for(int it = 0; it < MAXN; it++){
		diff = bd = 0;
		for(int i = 0; i < n; i++) {
			double las = x[i];
			x[i] = 0;
			for(int j = 0; j < n; j++) if(i != j) {
				x[i] -= a[i][j] * x[j];
			}
			x[i] += b[i];
			x[i] /= a[i][i];
			if(fabs(x[i] - las) > diff) diff = fabs(x[i] - las);
			if(fabs(x[i]) > bd) bd = fabs(x[i]);
		}
		if(bd > bound) return -2;
		if(diff < TOL) return it+1;
	}
	return 0;
}

int main()
{
	int n, MAXN, i, j, k;
	double a[MAX_SIZE][MAX_SIZE], b[MAX_SIZE], x[MAX_SIZE];
	double TOL;

	scanf("%d", &n);
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++)
			scanf("%lf", &a[i][j]);
		scanf("%lf", &b[i]);
	}
	scanf("%lf %d", &TOL, &MAXN);

	printf("Result of Jacobi method:\n");
	for ( i=0; i<n; i++ )
		x[i] = 0.0;
	k = Jacobi( n, a, b, x, TOL, MAXN );
	switch ( k ) {
		case -2:
			printf("No convergence.\n");
			break;
		case -1: 
			printf("Matrix has a zero column.  No unique solution exists.\n");
			break;
		case 0:
			printf("Maximum number of iterations exceeded.\n");
			break;
		default:
			printf("no_iteration = %d\n", k);
			for ( j=0; j<n; j++ )
				printf("%.8f\n", x[j]);
			break;
	}
	printf("Result of Gauss-Seidel method:\n");
	for ( i=0; i<n; i++ )
		x[i] = 0.0;
	k = Gauss_Seidel( n, a, b, x, TOL, MAXN );
	switch ( k ) {
		case -2:
			printf("No convergence.\n");
			break;
		case -1: 
			printf("Matrix has a zero column.  No unique solution exists.\n");
			break;
		case 0:
			printf("Maximum number of iterations exceeded.\n");
			break;
		default:
			printf("no_iteration = %d\n", k);
			for ( j=0; j<n; j++ )
				printf("%.8f\n", x[j]);
			break;
	}

	return 0;
}