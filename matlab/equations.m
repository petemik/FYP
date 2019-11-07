function output = equations(t, y)
output = zeros(4, 1);

xx = y(1);
yy = y(2);

m_sun = 1.989e30;
G = 6.674e-11;
% m_earth = 5.972e24;
r = (xx.^2+yy.^2).^0.5;

%year = np.sqrt((4*(np.pi^2)*a^3)/(G*m_sun));


output(1) = y(3);
output(3) = -G.*m_sun.*xx.*r.^-3;
output(2) = y(4);
output(4) = -G*m_sun.*yy.*r.^-3;
end

