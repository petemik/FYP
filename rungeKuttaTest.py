import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import RK45, Radau, BDF


# mass in earth mass
m_sun = 1.989e30
m_earth = 5.972e24
year = 3.154e7
#G = 6.674e-11
# In au^3 EM^-1 s^-2
G = 6.674e-11

## ISSUE: The issue appears to be due to the numbers getting too small at a certain point
def func(t, y):
    dxdt = y[2]
    dydt = y[3]
    dvxdt = y[0]*(-G*m_sun)/(np.power((np.power(y[0],2)+np.power(y[1],2)),3/2))
    dvydt = y[1]*(-G*m_sun)/(np.power((np.power(y[0],2)+np.power(y[1],2)),3/2))
    #dvydt = 0
    return [dxdt, dydt, dvxdt, dvydt]

answer = np.array([])
time = np.array([])
t0 = 0

# IN au
x0 = 150e9
y0 = 0 
vx0 = 0
vy0 = 30000

initial_state = [x0, y0, vx0, vy0]

x_vals = []
y_vals = []

max_steps=10000000
sol =  RK45(func, t0=t0, y0=initial_state, t_bound=10*year)
for i in range(max_steps):
    
    sol.step()
    print(sol.status)
    if sol.status=="failed" or sol.status=="finished":
        print("Integrator has failed on step", i)
        print(sol.y)
        break
    else: 
        time = np.append(time, sol.t)
        x_vals = np.append(x_vals, sol.y[0])
        y_vals = np.append(y_vals, sol.y[1])
        
#plt.plot(time, x_vals)
plt.plot(x_vals, y_vals)
#plt.plot(time, answer)