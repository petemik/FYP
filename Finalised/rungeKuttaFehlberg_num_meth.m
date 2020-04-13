function [time, sol5] =rungeKuttaFehlberg_num_meth(func, tspan, initial, max_stepsize, atol, rtol)
    %Initialisation
    j=1;
    sol4(:, 1) = initial';
    sol5(:, 1) = initial';
    time(1) = tspan(1);
    stepsize = max_stepsize;
    n_eqs = size(initial,1);
    S = 0.95;
  
    %Iterate through RKF method
    while time(j)<=tspan(2)
        err = 2;
        while (err > 1)
            k1 = stepsize*func(time(j)...
                , sol5(:, j));
            k2 = stepsize*func(time(j) + (1/5)*stepsize...
                , sol5(:,j)+(1/5).*k1);
            k3 = stepsize*func(time(j) + (3/10)*stepsize...
                , sol5(:, j)+(3/40).*k1+(9/40).*k2);
            k4 = stepsize*func(time(j) + (4/5)*stepsize...
                , sol5(:, j)+(44/45).*k1+(-56/15).*k2+(32/9).*k3);
            k5 = stepsize*func(time(j) + (8/9)*stepsize...
                , sol5(:, j)+(19372/6561).*k1+(-25360/2187).*k2+(64448/6561).*k3+(-212/729).*k4);
            k6 = stepsize*func(time(j) + stepsize...
                , sol5(:, j)+(9017/3168).*k1+(-355/33).*k2+(46732/5247).*k3+(49/176).*k4+(-5103/18656).*k5);
            sol4(:, j+1) = sol5(:, j) ...
                               +(5179/57600).*k1+(7571/16695).*k3+(393/640).*k4+(-92097/339200).*k5+(187/2100).*k6;
            sol5(:, j+1) = sol5(:, j) ...
                               +(35/384).*k1+(500/1113).*k3+(125/192).*k4+(-2187/6784).*k5+(11/84).*k6;
            
            
            scale = atol + rtol*max(norm(sol5(:,j+1)),norm(sol5(:,j)));
            triangle = abs(sol5(:,j+1)-sol4(:,j+1));
            triangleDivScale = triangle/scale;
            err = (n_eqs)^(-1/2)*norm(triangleDivScale);
            if err ~= 0
                stepsize = S*stepsize*(abs(err))^(-1/5);
            end
        end
        time(j+1) = time(j) + stepsize;
        j=j+1;
    end
        
end