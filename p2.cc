#include <stdio.h>
#include <math.h>

#define ZERO 1e-13 /* X is considered to be 0 if |X|<ZERO */
#define MAXN 11    /* Max Polynomial Degree + 1 */
double calc(int n, double f[], double a) {
    double ret = 0, x = 1;
    for(int i = 0; i <= n; i++, x *= a) {
        ret += f[i] * x;
    }
    return ret;
}
double getRoot(double start, int n, double f[], double df[]) { 
    double p = start, lp = 0;
    int cnt = 0;
    double num = 0, dom = 0, vf = 1e9, vdf, vddf;
    while(fabs(p - lp) > 1e-10) {
        cnt++;
        if(cnt > 10) break;
        lp = p;
        num = calc(n, f, lp);
        dom = calc(n-1, df, lp);
        if (fabs(dom) < 1e-12) {
            break;
        }
        p = lp - num / dom;
        // printf("%lf  %lf    %lf\n",num, dom, p);
    }
    return p;
}
double getRoot2(double start, int n, double f[], double df[], double ddf[]) { 

    double p = start, lp = 0;
    int cnt = 0;
    double num = 0, dom = 0, vf = 1e9, vdf, vddf;
    while(fabs(p - lp) > 1e-10) {
        cnt++;
        if(cnt > 10) break;
        lp = p;
        vf = calc(n, f, lp);
        vdf = calc(n-1, df, lp);
        vddf = calc(n-2, ddf, lp);
        num = vf * vdf;
        dom = vdf*vdf - vf*vddf;
        if (fabs(dom) < 1e-12) {
            break;
        }
        p = lp - num / dom;
        // printf("%lf  %lf    %lf\n",num, dom, p);
    }
    return p;
}
double Polynomial_Root(int n, double f[], double a, double b, double EPS) {
    double df[MAXN] = {0};
    for(int i = 0; i < n; i++) {
        df[i] = f[i+1] * (i+1);
    }
    double ddf[MAXN] = {0};
    for(int i = 0; i < n-1; i++) {
        ddf[i] = df[i+1] * (i+1);
    }
    double ret = a, minv = 1e9, v, p;
    ret = getRoot2((a+b)/2, n, f, df, ddf);
    minv = calc(n, f, ret);
    p = getRoot((a+b)/2, n, f, df);
    v = calc(n, f, p);
    if(fabs(minv > v)) ret = p;
    return ret;
}

int main()
{
    int n;
    double c[MAXN], a, b;
    double EPS = 0.00005;
    int i;

    scanf("%d", &n);
    for (i=n; i>=0; i--) 
        scanf("%lf", &c[i]);
    scanf("%lf %lf", &a, &b);
    printf("%.4f\n", Polynomial_Root(n, c, a, b, EPS));

    return 0;
}