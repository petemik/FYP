t = csvread('RK_t_432000_321.6624.csv');
% t=t';
iterations = size(t,1);
H=zeros(iterations-1,1);

for i=1:iterations-1
    H(i,1) = t(i+1,1) - t(i,1);
end
plot(H)