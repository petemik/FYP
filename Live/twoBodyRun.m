% Define constants
year = 3.154e7;
millenia = 1000*year;

%planets initial conditions [x0;y0;vx0;vy0] & masses
sun =       [0;         0;  0;  0];
%earth =     [147.1e9;     0;  0;  30.3e3];
jupiter =   [740e9;    0;  0;  13.1e3];
masses =    [2e30; 
            1898e24];
initial = [ sun;
            jupiter];

tic
%Perform in-built Runge-Kutta
func = @(t, y) multiEqs(t, y, masses);
num_years = 1000;
end_time = num_years*year;
tspan = [0 end_time];
stepsize = (1/1000)*year;
% [t,y] = rungeKuttaFehlberg56(func, tspan, initial, stepsize, 1e-6);
% [t,y] = rungeKutta(func, tspan, initial, stepsize);
opts = odeset('RelTol',1e-2,'AbsTol',1e-4, 'MaxStep', stepsize);
[t,y] = ode45(func,tspan, initial, opts);
y = transpose(y);
toc

%Plot orbits 
steps = max(size(t));
% Earth Orbit analytical

rmax = 152.1e9;
rmin = 147.1e9;
E = (rmax-rmin)/(rmax+rmin);
a = 149.60e9; % Semi-major axis
bsquared = (rmax*rmin); % Semi-minor axis
p = bsquared/a;
theta = 0:2*num_years*pi/steps:(2*num_years*pi-2*num_years*pi/steps);
steps = size(theta, 2);
r = zeros(1, steps);
for i=1:steps
    r(i) = a*(1-E^2)/(1+E*cos(theta(i)));
end
x_anal = r.*cos(theta);
y_anal = r.*sin(theta);

r_model = (y(4*1+1,:).^2 + y(4*1+2, :).^2).^(1/2);
diff_r = r - r_model;
hold on
n_bodies = size(masses, 1);
plot(y(1,:), y(2, :), '-')
plot(y(4*1+1,:), y(4*1+2, :), 'r', x_anal, y_anal,'b');

plot(diff_r);
% legend('Sun', 'Earth (Model)', 'Earth (Analytical)')
hold off




