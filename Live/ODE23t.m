% Define constants
G = 6.674e-11;
year = 3.154e7;
millenia = 1000*year;

%planets initial conditions [x0;y0;vx0;vy0] & masses
sun =       [0;         0;  0;  0];
mercury =   [46e9;      0;  0;  47.4e3];
venus =     [108e9;     0;  0;  35.3e3];
earth =     [147e9;     0;  0;  30.3e3];
mars =      [206e9;     0;  0;  26.5e3];
jupiter =   [740e9;     0;  0;  13.1e3];
saturn =    [135e10;    0;  0;  9.7e3];
uranus =    [2.7e12;    0;  0;  6.8e3];
neptune =   [4.5e12;    0;  0;  5.4e3];
masses = [  2e30; 
            0.33e24;
            4.87e24;
            5.97e24;
            0.642e24;
            1898e24;
            568e24;
            86.7e24;
            102e24];
initial = [ sun;
            mercury; 
            venus; 
            earth; 
            mars; 
            jupiter; 
            saturn; 
            uranus; 
            neptune];

%Perform in-built Runge-Kutta
func = @(t, y) multiEqs(t, y, masses);
end_time = millenia;
tspan = [0 end_time];
opts = odeset('RelTol',1e-2,'AbsTol',1e-4, 'MaxStep', (1/1000)*year);
[t,y] = ode23t(func,tspan, initial, opts);
y = y.';
total_steps = size(y, 2);

n_bodies = size(masses, 1);
T = zeros(n_bodies, total_steps);
% Calculate the Kinetic Energy of each planet
for i=0:n_bodies-1
    T(i+1, :) = 1/2*masses(i+1)*(y(4*i+3, :).^2 + y(4*i+4, :).^2);
end
% Sum to a total kinetic energy of the system
total_T = sum(T);

L = zeros(n_bodies, total_steps);
% Calculate the angular momentum of each planet
for i=0:n_bodies-1
    L(i+1, :) = masses(i+1)*(y(4*i+3, :).*y(4*i+2, :) - y(4*i+4, :).*y(4*i+1, :));
end
% Sum to a total angular momentum of the system
total_L = sum(L);

% Calculate potential energy, this is gonna be long
% Loop round this and do it for every element of the sytem maybe who tf
% knows
V = zeros(n_bodies, n_bodies, total_steps);

% First lets just do potential energy of the Sun
for i=1:n_bodies
    for j=(i+1):n_bodies
        V(i, j, :) = -G*masses(i)*masses(j)*...
            ((y(4*(i-1)+1,:)- y(4*(j-1)+1, :)).^2 +(y(4*(i-1)+2, j, :)- y(4*(j-1)+2, :)).^2).^(-1/2);
    end
end

total_V = sum(V, [1, 2]);
total_V = reshape(total_V, [1, total_steps]);
E = (total_T + total_V);
delta_E = (E-E(1))./E(1);
% plot(delta_E*100);
delta_L = (total_L-total_L(1))./total_L(1);
plot(delta_L*100);
