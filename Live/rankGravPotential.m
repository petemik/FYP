masses =    [0.33e24;
             4.87e24;
             5.97e24;
             0.642e24;
             1898e24;
             568e24;
             86.7e24;
             102e24];

G = 6.674e-11;
MSun = 2e30;

mean_distance = [57e9;
                 108e9;
                 150e9;
                 228e9;
                 779e9;
                 1430e9;
                 2880e9;
                 4500e9];
num_planets = size(masses, 1);
grav_impact = zeros(1, num_planets);
force_impact = zeros(1, num_planets);
for i=1:num_planets
    grav_impact(i) = -G*masses(i)*MSun/(mean_distance(i));
    force_impact(i) = G*masses(i)*MSun/((mean_distance(i))^2);
end

             