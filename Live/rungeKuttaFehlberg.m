function [time, sol4] =rungeKuttaFehlberg(func, tspan, initial, max_stepsize, tol)
    stepsize = max_stepsize;
    n_eqs = size(initial, 1);
    total_steps = ceil((tspan(2)-tspan(1))/stepsize);
    time = zeros(1, total_steps);
    sol4 = zeros(n_eqs, total_steps);
    sol5 = zeros(n_eqs, total_steps);
    sol4(:, 1) = initial';
    sol5(:, 1) = initial';
    time(1) = tspan(1);
    j=1;
    current_time=tspan(1);
    
    while current_time<=tspan(2)
        k1 = stepsize*func(time(j)...
            , sol4(:, j));
        k2 = stepsize*func(time(j) + (1/4)*stepsize...
            , sol4(:,j)+(1/4).*k1);
        k3 = stepsize*func(time(j) + (3/8)*stepsize...
            , sol4(:, j)+(3/32).*k1+(9/32).*k2);
        k4 = stepsize*func(time(j) + (12/13)*stepsize...
            , sol4(:, j)+(1932/2197).*k1+(-7200/2197).*k2+(7296/2197).*k3);
        k5 = stepsize*func(time(j) + stepsize...
            , sol4(:, j)+(439/216).*k1+(-8).*k2+(3680/513).*k3+(-845/4104).*k4);
        k6 = stepsize*func(time(j) + (1/2)*stepsize...
            , sol4(:, j)+(-8/27).*k1+(2).*k2+(-3544/2565).*k3+(1859/4104).*k4+(-11/40).*k5);
        sol4(:, j+1) = sol4(:, j) ...
                           +(25/216).*k1+(1408/2565).*k3+(2197/4101).*k4+(-1/5).*k5;
        sol5(:, j+1) = sol4(:, j) ...
                           +(16/135).*k1+(6656/12825).*k3+(28561/56430).*k4+(-9/50).*k5+(2/55).*k6;
        time(j+1) = time(j) + stepsize;
        
        s = (tol*stepsize/(2*norm(sol5(:, j+1)-sol4(:, j+1)))).^(1/4);
        stepsize = min([s*stepsize max_stepsize]);
        current_time=time(j);
        j=j+1;
    end
        
end