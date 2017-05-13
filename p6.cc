#include <stdio.h>

#define MAX_N 10

void Thomas(int n, double A[][3], double x[], double b[]) {
	double L[n+2];
	double U[n+2][2];
	double y[n+2];
	for(int i = 0; i < n+2; i++) {
		L[i] = U[i][0] = U[i][1] = y[i] = 0;
	}
	U[0][0] = A[0][1];
	U[0][1] = A[0][2];
	for(int i = 1; i < n; i++) {
		L[i] = A[i][0] / U[i-1][0];
		U[i][0] = A[i][1] - A[i-1][2] * L[i];
		U[i][1] = A[i][2];
	}
	U[n-1][1] = 0;
	y[0] = b[0];
	for(int i = 1; i < n; i++) {
		y[i] = b[i] - y[i-1] * L[i];
	}
	x[n-1] = y[n-1] / U[n-1][0];
	for(int i = n-2; i >= 0; i--) {
		x[i] = (y[i] - U[i][1] * x[i+1]) / U[i][0];
	}
}
double dif(double f[], double x[], int i, int j) {
	return (f[i] - f[j]) / (x[i] - x[j]);
}
void Cubic_Spline(int n, double x[], double f[], int Type, double s0, double sn, double a[], double b[], double c[], double d[]) {
	double h[n+2];
	double g[n+2];
	double M[n+2];
	double A[n+2][3];
	for(int i = 0; i <= n; i++) {
		h[i] = g[i] = M[i] = 0;
		A[i][0] = A[i][1] = A[i][2] = 0;
	}
	for(int i = 1; i <= n; i++)
		h[i] = x[i] - x[i-1];
	for(int i = 1; i < n; i++) {
		A[i][2] = h[i+1] / (h[i] + h[i+1]);
		A[i][1] = 2;
		A[i][0] = 1 - A[i][2];
		g[i] = (dif(f, x, i, i+1) - dif(f, x, i-1, i)) * 6 / (h[i] + h[i+1]);
	}
	if(Type == 1) {
		A[0][1] = 2;
		A[0][2] = 1;
		g[0] = (dif(f, x, 0, 1) - s0) * 6 / h[1];

		A[n][0] = 1;
		A[n][1] = 2;
		g[n] = (sn - dif(f, x, n-1, n)) * 6 / h[n];
	}
	if(Type == 2) {
		A[0][1] = 2;
		A[0][2] = 0;
		g[0] = 2 * s0;

		A[n][0] = 0;
		A[n][1] = 2;
		g[n] = 2 * sn;
	}
	// printf("mu\t");
	// for(int i = 0; i <= n; i++)
	// 	printf("%lf%c", A[i][0], "\t\n"[i == n]);
	// printf("2\t");
	// for(int i = 0; i <= n; i++)
	// 	printf("%lf%c", A[i][1], "\t\n"[i == n]);
	// printf("lmd\t");
	// for(int i = 0; i <= n; i++)
	// 	printf("%lf%c", A[i][2], "\t\n"[i == n]);
	// printf("h\t");
	// for(int i = 0; i <= n; i++)
	// 	printf("%lf%c", h[i], "\t\n"[i == n]);
	// printf("g\t");
	// for(int i = 0; i <= n; i++)
	// 	printf("%lf%c", g[i], "\t\n"[i == n]);
	
	Thomas(n+1, A, M, g);

	// printf("M\t");
	// for(int i = 0; i <= n; i++)
	// 	printf("%lf%c", M[i], "\t\n"[i == n]);
	
	for(int i = 1; i <= n; i++) {
		d[i] = (M[i] - M[i-1]) / (6 * h[i]);
		c[i] = M[i-1] / 2;
		double Aj = (f[i] - f[i-1]) / h[i] - (M[i] - M[i-1]) * h[i] / 6;
		b[i] = Aj - M[i-1] * h[i] / 2;
		a[i] = f[i-1];
		// printf("%lf\n", a[i]);
	}
}

double S( double t, double Fmax, int n, double x[], double a[], double b[], double c[], double d[] ) {
	for(int i = 1; i <= n; i++) {
		if(x[i-1] <= t && t <= x[i]) {
			double p = t - x[i-1];
			return a[i] + b[i] * p + c[i] * p * p + d[i] * p * p * p;
		}
	}
	return Fmax;
}

int main()
{
	int n, Type, m, i;
	double x[MAX_N], f[MAX_N], a[MAX_N], b[MAX_N], c[MAX_N], d[MAX_N];
	double s0, sn, Fmax, t0, tm, h, t;

	scanf("%d", &n);
	for (i=0; i<=n; i++) 
		scanf("%lf", &x[i]);
	for (i=0; i<=n; i++) 
		scanf("%lf", &f[i]);
	scanf("%d %lf %lf %lf", &Type, &s0, &sn, &Fmax);

	Cubic_Spline(n, x, f, Type, s0, sn, a, b, c, d);
	for (i=1; i<=n; i++)
		printf("%12.8e %12.8e %12.8e %12.8e \n", a[i], b[i], c[i], d[i]);

	scanf("%lf %lf %d", &t0, &tm, &m);
	h = (tm-t0)/(double)m;
	for (i=0; i<=m; i++) {
		t = t0+h*(double)i;
		printf("f(%12.8e) = %12.8e\n", t, S(t, Fmax, n, x, a, b, c, d));
	}

	return 0;
}
