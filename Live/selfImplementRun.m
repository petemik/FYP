% Define constants
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
end_time = 10*year;
tspan = [0 end_time];
stepsize = (1/1000)*year;
[t,y] = rungeKuttaFehlberg56(func, tspan, initial, stepsize,1e-6);
%[t,y] = rungeKutta(func, tspan, initial, stepsize);
toc

% fileID = fopen('test-data.txt','w');
% fprintf(fileID,'%6s|%6s|%6s\n','t','x','y');
% fprintf(fileID,'%12.8f', y);
% fclose(fileID);

%Plot orbits 
hold on
n_bodies = size(masses, 1);
for i=0:n_bodies-1
    plot(y(4*i+1,:), y(4*i+2, :))
end
legend('Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune')
hold off
