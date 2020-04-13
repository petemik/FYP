% Note this code is written assuming Polar coordinates 
clear;
% Define constants
year = 3.154e7;
millenia = 1000*year;
day = 60*60*24;
G = 6.67e-11;

M = 2e30;                                                   % Sun mass
mu = 5.97e24;                                               % Earth reduced mass

%planets initial conditions [x0;y0;vx0;vy0] & masses
sun =       [0;	0; 0; 0];
earth =     [147e9;	0; 0; 30.3e3];
masses =    [M; mu];
initial =   [sun; earth];

% [r_initial, theta_inital, r_dot_initial, theta_dot_initial]
R_initial = [0 0 0 0];                                      % Sun IC
r_initial = [earth(1) 0 0 earth(4)/earth(1)];               % Earth IC
theta_0 = r_initial(2) - R_initial(2);                      % Theta IC

%Values
L = mu*r_initial(1)^2 * r_initial(4);
k = -G*M*mu;
E = -G*M*mu/r_initial(1) ...
    + 1/2*mu*(r_initial(3)^2 + r_initial(1)^2 * r_initial(4)^2)...
    + 1/2*M*(R_initial(3)^2 + R_initial(1)^2 * R_initial(4)^2);
e = (1 + 2*E*L^2/(k^2*mu))^(1/2);

%Perform in-built Runge-Kutta
func = @(t, y) multiEqs(t, y, masses);
end_time = year;
tspan = [0 end_time];
stepsize = 40*day;  
H = stepsize;
steps = ceil((tspan(2)-tspan(1))/H);
kmax = 8;
t0 = tspan(1);
y(1, :) = initial;

tic
%     [t, y] = rungeKutta(func, tspan, initial, stepsize);
    [t, y] = BS1DAdaptivev2_2(func, tspan, initial, stepsize, kmax,1E-14,1E-12);
toc

% tic
% for i=1:steps
%     small_tspan = [t0+(i-1)*H t0+i*H];
%     y(i+1, :) = BS1DAdaptive(func, small_tspan, y(i, :), H, kmax, 1e-6, 1e-4);
% end
% duration = toc;
% y=y';
% t = 0:H:tspan(2);


%Numerical sol (r,theta)
y_numerical = zeros(2,size(y,2));
for i=1:size(y,2)
    y_numerical(1,i) = ((y(5,i)-y(1,i))^(2)+(y(6,i)-y(2,i))^(2))^(1/2);
    if (y(5,i) - y(1,i) > 0) && (y(6,i) - y(2,i) > 0) 
        y_numerical(2,i) = atan((y(6,i)-y(2,i))/(y(5,i)-y(1,i)));
    elseif (y(5,i) - y(1,i) > 0) && (y(6,i) - y(2,i) < 0)
        y_numerical(2,i) = 2*pi + atan((y(6,i)-y(2,i))/(y(5,i)-y(1,i)));
    elseif (y(5,i) - y(1,i) < 0) && (y(6,i) - y(2,i) > 0)
        y_numerical(2,i) = pi + atan((y(6,i)-y(2,i))/(y(5,i)-y(1,i)));
    elseif (y(5,i) - y(1,i) < 0) && (y(6,i) - y(2,i) < 0)
        y_numerical(2,i) = pi + atan((y(6,i)-y(2,i))/(y(5,i)-y(1,i)));
    end
end  

theta = y_numerical(2,:);

%Analytic soution
u = -k*mu/L^2 * (1 + e*cos(theta - theta_0));
r = 1./u;

%Error
y_error = abs((y_numerical(1,:)-r)./r);
plot(y_error);

max_error = max(y_error);




