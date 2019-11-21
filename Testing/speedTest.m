func = @(t, x) 3*t^2;

i=1;
initial = 0;
tspan = [0 1000];
stepsize=(1/100);

tic
    [t,y] = rungeKutta(func, tspan, initial, stepsize);
toc

A=transpose(t);
B=transpose(y);

fileID = fopen('speedtest-data.txt','w');
fprintf(fileID,'%6s|%6s\n','t','x');
for i=1:1000:1000001
    fprintf(fileID,'%12.8f|%12.8f\n', A(i),B(i));
end
fclose(fileID);