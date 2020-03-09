% Define constants
year = 3.154e7;
millenia = 1000*year;

%planets initial conditions [x0;y0;vx0;vy0] & masses
sun =       [-6.4382245E+08;	1.0921744E+09;	-1.4241571E+01;	-4.4343761E+00];
mercury =   [-5.7995526E+10;	2.8370562E+09;	-1.1628598E+04;	-4.6603509E+04];
venus =     [-2.3107878E+10;	1.0625732E+11;	-3.4381530E+04;	-7.5130530E+03];
earth =     [-1.4132709E+11;	4.7898844E+10;	-9.8949287E+03;	-2.8387273E+04];
mars =      [-9.8203828E+10;	-2.0338091E+11;	2.2767800E+04;	-8.3582128E+03];
jupiter =   [1.4542722E+11;     -7.6379101E+11;	1.2675763E+04;	3.0667018E+03];
saturn =    [6.1148506E+11;     -1.3685160E+12;	8.2848034E+03;	3.9134473E+03];
uranus =    [2.4057104E+12;     1.7310000E+12;	-4.0273876E+03;	5.2104067E+03];
neptune =   [4.3799648E+12; 	-9.2323714E+11;	1.0847279E+03;	5.3513090E+03];
masses =    [2e30; 
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

tic
%Perform in-built Runge-Kutta
func = @(t, y) multiEqs(t, y, masses);
end_time = millenia;
tspan = [0 end_time];
stepsize = (1/1000)*year;
%[t,y] = rungeKuttaFehlberg56(func, tspan, initial, stepsize,1e-6);
%[t,y] = rungeKutta(func, tspan, initial, stepsize);
opts = odeset('RelTol',1e-2,'AbsTol',1e-4, 'MaxStep', stepsize);
[t,y] = ode45(func,tspan, initial, opts);
y = transpose(y);
toc

% fileID = fopen('test-data.txt','w');
% fprintf(fileID,'%6s|%6s|%6s\n','t','x','y');
% fprintf(fileID,'%12.8f', y);
% fclose(fileID);

%Plot orbits 
hold on
n_bodies = size(masses, 1);
for i=0:n_bodies-1
    plot(y(4*i+1,:), y(4*i+2,:))
end
legend('Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune')
hold off
