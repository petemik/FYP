% Define constants
clear
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

%Perform in-built Runge-Kutta
f = @(t, y) multiEqs(t, y, masses);
end_time = millenia;
tspan = [0 end_time];
kmax = 8;
range = [5; 10; 20; 40];
n_eqs = size(initial, 1);

for i=1:1
    stepsize = range(i)*day;
    tic
        [t, y] = BS1DAdaptivev2(f, tspan, initial, stepsize, kmax,1E-14,1E-12);
    duration = toc;

    % n = 1;
    % y_slim = y(: , 1:n:end);
    % t_slim = t(: , 1:n:end);
    writematrix(y,'BSa_'+string(stepsize)+'_'+string(duration)+'.csv');
    writematrix(t,'BSa_t_'+string(stepsize)+'_'+string(duration)+'.csv');
end

sound(sin(1:3000));
