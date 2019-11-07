func = @(t, y) equations(t, y);

m_sun = 1.989e30;
m_earth = 5.972e24;
a=384400e3;
e=0.0167;

G = 6.674e-11;
%year = 4*(pi^2)*a^(1/2)/(G*m_sun);
year = 3.154e7;
millenia = 1000*year;

x0 = 149.513e9;
y0 = 0;
vx0 = 0;
vy0 = 29.78e3;

initial = [x0; y0; vx0; vy0];
tspan = [0 100000*31.536e6];
%opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
opts = odeset('RelTol',1e-2,'AbsTol',1e-4, 'MaxStep', (1/1000)*year);
[t,y] = ode45(func,tspan, initial, opts);

plot(y(:,1), y(:,2))

