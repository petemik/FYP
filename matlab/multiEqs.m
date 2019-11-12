function output = multiEqs(t, y)
n_bodies = size(y, 1)/4;
output = zeros(4*n_bodies, 1);
% Here xx and yy represent distance of each planet to the Sun.
xx = zeros(n_bodies, 1);
yy = zeros(n_bodies, 1);
for i = 0:n_bodies-1
    xx(i+1) = y(4*i + 1);
    yy(i+1) = y(4*i + 2);
end
M = 1.989e30;
m = 5e24;
G = 6.674e-11;
r = (xx.^2+yy.^2).^0.5;

% Here lets do a matrix of the different combinations
% In the future for efficiency we will just fill in half and then take the
% negative of it in the transpose (don't know if that makes sense)
Fx = zeros(n_bodies, n_bodies);
Fy = zeros(n_bodies, n_bodies);
for i=1:n_bodies
    for j=1:n_bodies
        if i ~= j
            Fx(i, j) = -G.*m.*(xx(i)-xx(j)).*(((xx(i)-xx(j)).^2 ...
                                            + (yy(i)-yy(j)).^2).^(-3/2));
            Fy(i, j) = -G.*m.*(yy(i)-yy(j)).*(((xx(i)-xx(j)).^2 ...
                                            + (yy(i)-yy(j)).^2).^(-3/2));
        else
        end
    end 
end

totalFx = sum(Fx, 2);
totalFy = sum(Fy, 2);

for i = 0:n_bodies-1
    output(4*i+1) = y(4*i+3);    
    output(4*i+3) = -G.*M.*xx(i+1).*r(i+1).^-3 + totalFx(i+1);
    output(4*i+2) = y(4*i + 4);
    output(4*i+4) = -G*M.*yy(i+1).*r(i+1).^-3 + totalFy(i+1);
end
end

