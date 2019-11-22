#include <iostream>
#include <math.h>
#include<stdlib.h>
#include <vector> 

using namespace std;

// Above here we need the multiEqs function which will be passed in as f to the rungeKutta

double rungeKutta(double(*f)(double, double, double), double tspan[], double initial[], double stepsize) {
	int n_eqs = sizeof(initial) / sizeof(double);
	int total_steps = ceil((tspan[1] - tspan[0]) / stepsize);
	double *time = (double *)malloc(total_steps * sizeof(double));
	for (int i = 0; i < total_steps; i++) {
		printf("%d", time[i]);
	}
	memset(time, 0, sizeof(time));
	return 7.3;
}

double * multiEqs(double t, double y[], double masses[]) {
	int n_bodies = sizeof(y) / (4 * sizeof(double));
	double* output = new double[4 * n_bodies];
	double* xx = new double[n_bodies];
	double* yy = new double[n_bodies];
	for (int i = 0; i < n_bodies; i++) {
		xx[i] = y[4 * i];
		yy[i] = y[4 * i + 1];
	}

	double G = 6.674*pow(10, -11);

	double** Fx = new double*[n_bodies];
	double** Fy = new double*[n_bodies];
	for (int i = 0; i < n_bodies; i++) {
		Fx[i] = new double[n_bodies];
		Fy[i] = new double[n_bodies];
	}

	for (int i = 0; i < n_bodies; i++) {
		for (int j = i + 1; j < n_bodies; i++) {
			Fx[i][j] = -G * masses[i] * masses[j] * (xx[i] - xx[j])*pow(pow(xx[i] - xx[j], 2) + pow(yy[i] - yy[j], 2), -3.0 / 2.0);
			Fx[i][j] = -G * masses[i] * masses[j] * (yy[i] - yy[j])*pow(pow(xx[i] - xx[j], 2) + pow(yy[i] - yy[j], 2), -3.0 / 2.0);
		}
	}

	for (int i = 0; i < n_bodies; i++) {
		for (int j = i + 1; j < n_bodies; i++) {
			Fx[j][i] = -Fx[i][j];
			Fy[j][i] = -Fy[i][j];
		}
	}

}

int main() {
	// Constants
	double G = 6.674*pow(10, -11);
	double year = 3.154*pow(10, 7);
	double millenia = 1000 * year;


	// Planets initial conditions [x0, y0, vx0, vy0] and masses
	int n_bodies = 9;
	double sun[4] = { 0, 0, 0, 0 };
	double mercury[4] = { 46 * pow(10, 9), 0, 0, 47.4*pow(10, 3) };
	double venus[4] = { 108 * pow(10, 9), 0, 0, 35.3*pow(10, 3) };
	double earth[4] = { 147 * pow(10, 9), 0, 0, 30.3*pow(10, 3) };
	double mars[4] = { 206 * pow(10, 9), 0, 0, 26.5*pow(10, 3) };
	double jupiter[4] = { 740 * pow(10, 9), 0, 0, 13.1*pow(10, 3) };
	double saturn[4] = { 135 * pow(10, 10), 0, 0, 9.7*pow(10, 3) };
	double uranus[4] = { 2.7*pow(10, 12), 0, 0, 6.8*pow(10, 3) };
	double neptune[4] = { 4.5*pow(10, 12), 0, 0, 5.4*pow(10, 3) };

	double* masses = new double[n_bodies];

	masses[0] = 1.99 * pow(10, 30);
	masses[1] = 0.33 * pow(10, 24);
	masses[2] = 4.87 * pow(10, 24);
	masses[3] = 5.97 * pow(10, 24);
	masses[4] = 0.64 * pow(10, 24);
	masses[5] = 1898 * pow(10, 24);
	masses[6] = 568 * pow(10, 24);
	masses[7] = 86.7 * pow(10, 24);
	masses[8] = 102 * pow(10, 24);

	double* initial = new double[4 * n_bodies];
	for (int i = 0; i < 4; i++) {
		initial[i] = sun[i];
		initial[i + 4] = mercury[i];
		initial[i + 8] = venus[i];
		initial[i + 12] = earth[i];
		initial[i + 16] = mars[i];
		initial[i + 20] = jupiter[i];
		initial[i + 24] = saturn[i];
		initial[i + 28] = uranus[i];
		initial[i + 32] = neptune[i];
	}

	double end_time = 100 * year;
	double tspan[2] = { 0, end_time };



	int dummy;
	cin >> dummy;
	return 0;

}