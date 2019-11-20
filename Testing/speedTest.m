func = @(t, x) 3*t^2;

initial = 0;
tspan = [0 1000];
opts = odeset('RelTol',1e-2,'AbsTol',1e-4, 'MaxStep', 1e-4);
[t,y] = ode45(func,tspan, initial, opts);
plot(t,y)
% plot(t,t.^3 - y)