// C program to implement Runge Kutta method 
#include<stdio.h> 
#include<math.h> 

float func(float x[4], float t) 
{ 
	 float m_sun = 1.989e30;
	 float G = 6.674e-11;
	 
	 float dxdt=x[2];
	 float dydt=x[3];
	 float dvxdt=(x[0]*(-1)*G*m_sun)/(pow(pow(x[0],2)+pow(x[1],2)),(3/2));
	 float dvydt=(x[1]*(-1)*G*m_sun)/(pow(pow(x[0],2)+pow(x[1],2)),(3/2));
	 
	 return [dxdt,dydt,dvxdt,dvydt];
} 

float rungeKutta(float initial[4], float t0, int n) 
{ 

	int h=157700;
	float k1[4], k2[4], k3[4], k4[4]; 

	// Iterate for number of iterations 
	float solution[4] = initial; 
	float t=t0;
	for (int i=1; i<=n; i++) 
	{ 
		// Apply Runge Kutta Formulas to find 
		// next value of solution
		for(int j=0;j<4;j++)
		{
			k1[j] = h*func(solution[j], t)[j]; 
			k2[j] = h*func(solution[j] + 0.5*k1[j],t + 0.5*h)[j]; 
			k3[j] = h*func(solution[j] + 0.5*k2[j],t + 0.5*h)[j]; 
			k4[j] = h*func(solution[j] + k3[j],t + h)[j]; 
		}
		
		// Update next value of solution 
		for(int j=0;j<4;j++)
		{
			solution[j] = solution[j] + (1.0/6.0)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
		}

		// Update next value of t 
		t = t + h; 
	} 

	return solution; 
} 

// Driver method 
int main() 
{ 
	float initial=[150e9,0,0,30000], t0=0; 
	int n=200;
	
	printf("\nThe value of Earth state is : %f\n", 
			rungeKutta(initial,t0,n)); 
	return 0; 
} 
