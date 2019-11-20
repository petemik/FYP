function [time, solution] =rungeKutta(func, tspan, initial, stepsize)
    n_eqs = size(initial, 1);
    total_steps = ceil((tspan(2)-tspan(1))/stepsize);
    time = zeros(1, total_steps);
    solution = zeros(n_eqs, total_steps);
    solution(:, 1) = initial';
    time(1) = tspan(1);
    k1=zeros(n_eqs*4,1);
    k2=zeros(n_eqs*4,1);
    k3=zeros(n_eqs*4,1);
    k4=zeros(n_eqs*4,1);
    
        for j=1:total_steps
            k1 = stepsize*func(time(j), solution(:, j));
            k2 = stepsize*func(time(j) + 0.5*stepsize, solution(:,j)+0.5.*k1);
            k3 = stepsize*func(time(j) + 0.5*stepsize, solution(:, j)+0.5.*k2);
            k4 = stepsize*func(time(j) + stepsize, solution(:, j) + k3);
            solution(:, j+1) = solution(:, j) ...
                               + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
            time(j+1) = time(j) + stepsize;
        end
        
end