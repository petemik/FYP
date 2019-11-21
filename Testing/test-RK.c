#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double rk4(double(*f)(double, double), double dt, double t, double x)
{
	double	k1 = dt * f(t, x),
		k2 = dt * f(t + dt / 2, x + k1 / 2),
		k3 = dt * f(t + dt / 2, x + k2 / 2),
		k4 = dt * f(t + dt, x + k3);
	return x + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

double rate(double t, double x)
{
	return 3*pow(t,2);
}

int main(void)
{
	clock_t start, end;
	double cpu_time_used;

	double* x, t, x2;
	double t0 = 0, t1 = 1000, dt = .01;
	int i, n = 1 + (t1 - t0) / dt;
	x = (double*)malloc(sizeof(double) * n);

	start = clock();
	for (x[0] = 0, i = 1; i < n; i++)
		x[i] = rk4(rate, dt, t0 + dt * (i - 1), x[i - 1]);
	end = clock();

	printf("t\tx.\n------------\n");
	for (i = 0; i < n; i += 1000) {
		t = t0 + dt * i;
		x2 = pow(t * t / 4 + 1, 2);
		printf("%lf\t%lf\n", t, x[i]);
	}

	cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
	printf("Runge-Kutta time taken: %lf\n", cpu_time_used);

	return 0;
}