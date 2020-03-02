% Define constants
year = 3.154e7;
millenia = 1000*year;

%planets initial conditions [x0;y0;vx0;vy0] & masses
sun =       [0;         0;  0;  0];

jupiter =   [740e9;    0;  0;  13.72e3];
masses =    [2e30; 
            1898e24];
initial = [ sun;
            jupiter];

tic
%Perform in-built Runge-Kutta
func = @(t, y) multiEqs(t, y, masses);
end_time = millenia;
tspan = [0 end_time];
stepsize = (1/1000)*year;
% [t,y] = rungeKuttaFehlberg56(func, tspan, initial, stepsize,1e-6);
opts = odeset('RelTol',1e-2,'AbsTol',1e-4, 'MaxStep', stepsize);
[t,y] = ode45(func,tspan, initial, opts);
y = transpose(y);
%[t,y] = rungeKutta(func, tspan, initial, stepsize);
toc

% fileID = fopen('test-data.txt','w');
% fprintf(fileID,'%6s|%6s|%6s\n','t','x','y');
% fprintf(fileID,'%12.8f', y);
% fclose(fileID);

%Plot orbits 
hold on
n_bodies = size(masses, 1);
% for i=0:n_bodies-1
%     plot(y(4*i+1,:), y(4*i+2, :))
% end
rSun = (y(4*0+1,:).^2 + y(4*0+2, :).^2).^(1/2);
rJup = ((y(4*1+1,:)-y(4*0+1,:)).^2 + (y(4*1+2, :)-y(4*0+2, :)).^2).^(1/2);
plot(y(4*1+1,1:400065), y(4*1+2, 1:400065))
legend('Sun', 'Jupiter')
hold off
