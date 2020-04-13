% Define constants
clear;
G = 6.674e-11;
year = 3.154e7;
millenia = 1000*year;
day = 60*60*24;

%planets initial conditions [x0;y0;vx0;vy0] & masses (as of 01-01-2020)
sun =       [-5.6827347E+08;	1.1130130E+09;      -1.4461742E+01;     -3.4475562E+00];
mercury =   [-1.0043171E+10;	-6.7829448E+10;     3.8473199E+04;      -4.1587487E+03];
venus =     [1.0762249E+11; 	8.9742506E+09;      -2.6935234E+03;     3.4766999E+04];
earth =     [-2.5453599E+10;	1.4609342E+11;      -2.9863807E+04;     -5.1658958E+03];
mars =      [-1.9805637E+11;	-1.3139630E+11;     1.4392944E+04;      -1.8050298E+04];
jupiter =   [7.8143326E+10;     -7.7695995E+11;     1.2839493E+04;      1.9303845E+03];
saturn =    [5.6749956E+11; 	-1.3883862E+12;     8.4066132E+03;      3.6278265E+03];
uranus =    [2.4267661E+12; 	1.7034167E+12;      -3.9624080E+03; 	5.2565982E+03];
neptune =   [4.3741565E+12; 	-9.5141866E+11;     1.1188384E+03;      5.3427212E+03];
masses =    [2e30; 0.33e24;     4.87e24; 5.97e24; 0.642e24; 1898e24; 568e24; 86.7e24; 102e24];
initial = [ sun;
            mercury; 
            venus; 
            earth; 
            mars; 
            jupiter; 
            saturn; 
            uranus; 
            neptune];

%Perform numerical method
% func = @(t, y) multiEqs(t, y, masses);
% end_time = 10*year;
% kmax = 8;
% tspan = [0 end_time];
% stepsize = 5*day;

% tic
% % [t, y] = rungeKuttaFehlberg_num_meth(func, tspan, initial, stepsize,1E-8,1E-6);
% [t, y] = BS1DAdaptivev2(func, tspan, initial, stepsize, kmax,1E-10,1E-8);
% % [t,y] = rungeKutta(func,tspan, initial, stepsize);
% % [t,y] = rungeKuttaFehlberg(func, tspan, initial, stepsize,1e-6);
% 
% % opts = odeset('RelTol',1e-2,'AbsTol',1e-4, 'MaxStep', stepsize);
% % [t,y] = ode45(func,tspan, initial, opts);
% % y = transpose(y);
% toc

y1 = csvread('BSa_t_432000_1901.6616.csv');
% y=y1';
% a = size(y1,1)/36;
% y_test=reshape(y1,[a,36]);
% y=y_test';
y=y1;

%data recording variables
increment = 1;
total_steps = size(y, 2);
recorded_steps = ceil(total_steps/increment);
n_bodies = size(masses, 1);

% Calculate the Kinetic Energy of each planet
T = zeros(n_bodies, recorded_steps);
for i=0:n_bodies-1
    for k=1:increment:total_steps
        index = ceil(k/increment);
        T(i+1, index) = 1/2*masses(i+1)*(y(4*i+3, k).^2 + y(4*i+4, k).^2);
    end
end
% Sum to a total kinetic energy of the system
total_T = sum(T);

% Calculate potential energy
V = zeros(n_bodies, n_bodies, recorded_steps);
for i=1:n_bodies
    for j=(i+1):n_bodies
        for k=1:increment:total_steps
            index = ceil(k/increment);
            V(i, j, index) = -G*masses(i)*masses(j)*...
                ((y(4*(i-1)+1,k)- y(4*(j-1)+1, k)).^2 +(y(4*(i-1)+2, k)- y(4*(j-1)+2, k)).^2).^(-1/2);
        end
    end
end
% Sum to a total potential energy of the system
total_V = sum(V, [1, 2]);
total_V = reshape(total_V, [1, recorded_steps]);

% Calculate the angular momentum of each planet
L = zeros(n_bodies, recorded_steps);
for i=0:n_bodies-1
    for k=1:increment:total_steps
        index = ceil(k/increment);
        L(i+1, index) = masses(i+1)*(y(4*i+3, k).*y(4*i+2, k) - y(4*i+4, k).*y(4*i+1, k));
    end
end
% Sum to a total angular momentum of the system
total_L = sum(L);

%Calculate percentage divergence from initial value
E = (total_T + total_V);
percentage_E = ((E-E(1))./E(1))*100;
percentage_L = ((total_L-total_L(1))./total_L(1))*100;

%plots
% figure(1);
% plot(percentage_E);
% figure(2);
% plot(percentage_L);

analysis(1,1)=max(percentage_E);
analysis(2,1)=min(percentage_E);
analysis(3,1)=mean(percentage_E);
analysis(4,1)=range(percentage_E);
analysis(1,2)=max(percentage_L);
analysis(2,2)=min(percentage_L);
analysis(3,2)=mean(percentage_L);
analysis(4,2)=range(percentage_L);

