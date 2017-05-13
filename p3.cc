#include <bits/stdc++.h>
using namespace std;
#define Max_size 100 /* max number of dishes */
void Thomas(int n, double A[][3], double x[], double b[]) {
    double L[Max_size] = {0};
    double U[Max_size][2] = {0};
    double y[Max_size] = {0};
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
double dot(int n, double a[], double b[]) {
    double ret = 0;
    for(int i = 0; i < n; i ++)
        ret += a[i] * b[i];
    return ret;
}
void Price( int n, double p[] ) {
	double A[Max_size][3] = {0};
    double u[Max_size] = {0};
    double v[Max_size] = {0};
    double y[Max_size] = {0};
    double q[Max_size] = {0};
	for(int i = 0; i < n; i++) {
		A[i][0] = 0.5;
        A[i][1] = 2;
		A[i][2] = 0.5;
	}
    u[0] = -A[0][1];
    u[n-1] = A[n-1][2];
    v[0] = 1;
    v[n-1] = -A[0][0] / A[0][1];
    A[n-1][1] -= u[n-1] * v[n-1];
    A[n-1][2] -= u[n-1] * v[0];
    A[0][1] -= u[0] * v[0];
    A[0][0] -= u[0] * v[n-1];
    Thomas(n, A, y, p);
    Thomas(n, A, q, u);
    double k = (dot(n, v, y) / (1 + dot(n, v, q)));
    for(int i = 0; i < n; i++) p[i] = y[i] - k * q[i];
    for(int i = 0; i < n; i++) printf("%.2lf ",y[i]); puts("");
    for(int i = 0; i < n; i++) printf("%.2lf ",q[i]); puts("");
}
int main()
{
    
    int n, i;
    double p[Max_size];

    scanf("%d", &n);
    for (i=0; i<n; i++) 
        scanf("%lf", &p[i]);
    Price(n, p);
    for (i=0; i<n; i++)
        printf("%.2f ", p[i]);
    printf("\n");

    return 0;
}