function output = multiEqs(t, y, masses)
n_bodies = size(y, 1)/4;
output = zeros(4*n_bodies, 1);
% Here xx and yy represent distance of each planet to the Sun.
xx = zeros(n_bodies, 1);
yy = zeros(n_bodies, 1);
for i = 0:n_bodies-1
    xx(i+1) = y(4*i + 1);
    yy(i+1) = y(4*i + 2);
end
G = 6.674e-11;

% Calculate forces
Fx = zeros(n_bodies, n_bodies);
Fy = zeros(n_bodies, n_bodies);
for i=1:n_bodies
    for j=(i+1):n_bodies
            Fx(i, j) = -G.*masses(i).*masses(j).*(xx(i)-xx(j)).*(((xx(i)-xx(j)).^2 ...
                                            + (yy(i)-yy(j)).^2).^(-3/2));
            Fy(i, j) = -G.*masses(i).*masses(j).*(yy(i)-yy(j)).*(((xx(i)-xx(j)).^2 ...
                                            + (yy(i)-yy(j)).^2).^(-3/2));
    end 
end

%Use transpose for reverse forces
Fx = Fx - transpose(Fx);
Fy = Fy - transpose(Fy);

totalFx = sum(Fx, 2);
totalFy = sum(Fy, 2);

for i = 0:n_bodies-1
    output(4*i+1) = y(4*i+3);    
    output(4*i+3) = totalFx(i+1)/masses(i+1);
    output(4*i+2) = y(4*i + 4);
    output(4*i+4) = totalFy(i+1)/masses(i+1);
end
end