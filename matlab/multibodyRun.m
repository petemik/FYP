% Define Constants
G = 6.674e-11;
year = 3.154e7;
millenia = 1000*year;

%planets initial conditions [x0;y0;vx0;vy0]
mercury =   [46e9; 0;  0;  47.4e3];
venus =     [108e9; 0;  0;  35.3e3];
earth =     [147e9; 0;  0;  30.3e3];
mars =      [206e9; 0;  0;  26.5e3];
jupiter =   [740e9; 0;  0;  13.1e3];
saturn =    [135e10;0;  0;  9.7e3];
uranus =    [2.7e12;0;  0;  6.8e3];
neptune =   [4.5e12;0;  0;  5.4e3];
masses = [  2e30; 
            0.33e24;
            4.87e24;
            5.97e24;
            0.642e24;
            1898e24;
            568e24;
            86.7e24;
            102e24];
initial = [ mercury; 
            venus; 
            earth; 
            mars; 
            jupiter; 
            saturn; 
            uranus; 
            neptune];

func = @(t, y) multiEqs(t, y, masses);
end_time = 100*year;
tspan = [0 end_time];
opts = odeset('RelTol',1e-2,'AbsTol',1e-4, 'MaxStep', (1/1000)*year);
[t,y] = ode45(func,tspan, initial, opts);
hold on
plot(0,0,'r*')
% Could do this in a loop
plot(y(:,1), y(:,2))
plot(y(:,5), y(:,6))
plot(y(:,9), y(:,10))
plot(y(:,13), y(:,14))
plot(y(:,17), y(:,18))
plot(y(:,21), y(:,22))
plot(y(:,25), y(:,26))
plot(y(:,29), y(:,30))
legend('Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune')
hold off
