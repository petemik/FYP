% Define constants
year = 3.154e7;
millenia = 1000*year;

%planets initial conditions [x0;y0;vx0;vy0] & masses (as of 01-01-2020)
wolf1061 =    [-5.146870672E+07;0;0;0];
wolf1061b =   [5.558531293E+09;0;0;8.33E+04];
wolf1061c =   [1.324853129E+10;0;0;5.41E+04];
wolf1061d =   [7.024853129E+10;0;0;2.35E+04];
masses =      [5.88E+29;	1.14027E+25;	2.03577E+25;	4.5969E+25];
initial = [ wolf1061;
            wolf1061b; 
            wolf1061c; 
            wolf1061d];

%Perform in-built Runge-Kutta
tic
func = @(t, y) multiEqs(t, y, masses);
end_time = 100*year;
tspan = [0 end_time];
stepsize = (1/10000)*year;
%[t,y] = rungeKuttaFehlberg56(func, tspan, initial, stepsize,1e-6);
% [t,y] = rungeKutta(func, tspan, initial, stepsize);
opts = odeset('RelTol',1e-2,'AbsTol',1e-4, 'MaxStep', stepsize);
[t,y] = ode45(func,tspan, initial, opts);
y = y';
toc

% Plot orbits 
hold on
n_bodies = size(initial, 1)/4;
for i=0:n_bodies-1
    plot(y(4*i+1,:), y(4*i+2,:))
end
legend('Wolf 1061', 'Wolf 1061b', 'Wolf 1061c', 'Wolf 1061d')
hold off

