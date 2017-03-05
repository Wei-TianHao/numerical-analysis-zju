#include <bits/stdc++.h>
using namespace std;
double l[100000];
void Series_Sum( double l[] ){
	for(int i = 0; i < 11; i += 1) {
		l[i] = 1;
		for(double j = 1; j < 200000; j += 1) {
			l[i] += (1-i/10.0)/(j*(j+1)*(j+(i/10.0)));
		}
	}
	for(int i = 11; i < 3001; i += 1) {
		double b = ((i-1) % 10 + 1) * 0.1;
		l[i] = l[((i-1) % 10 + 1)] * b;
		for(int j = 1; j <= (i-1) / 10; j ++)
			l[i] += 1.0/(j+b);
		l[i] = l[i] / ((i-1)/10 + b);
	}
}