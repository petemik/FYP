function [time, sol5] =rungeKuttaFehlberg56(func, tspan, initial, max_stepsize, tol)
   %Initialisation
    j=1;
    sol5(:, 1) = initial';
    sol6(:, 1) = initial';
    time(1) = tspan(1);
    stepsize = max_stepsize;
   
    %Iterate through RKF method
    while time(j)<=tspan(2)
        k1 = stepsize*func(time(j)...
            , sol5(:, j));
        k2 = stepsize*func(time(j) + (1/6)*stepsize...
            , sol5(:,j)+(1/6).*k1);
        k3 = stepsize*func(time(j) + (4/15)*stepsize...
            , sol5(:, j)+(4/75).*k1+(16/75).*k2);
        k4 = stepsize*func(time(j) + (2/3)*stepsize...
            , sol5(:, j)+(5/6).*k1+(-8/3).*k2+(5/2).*k3);
        k5 = stepsize*func(time(j) + (4/5)*stepsize...
            , sol5(:, j)+(-8/5).*k1+(144/25).*k2+(-4).*k3+(16/25).*k4);
        k6 = stepsize*func(time(j) + stepsize...
            , sol5(:, j)+(361/320).*k1+(-18/5).*k2+(407/128).*k3+(-11/80).*k4+(55/128).*k5);
        k7 = stepsize*func(time(j)...
            , sol5(:, j)+(-11/640).*k1+(11/256).*k3+(-11/160).*k4+(11/256).*k5);
        k8 = stepsize*func(time(j) + stepsize...
            , sol5(:, j)+(93/640).*k1+(-18/5).*k2+(803/256).*k3+(-11/160).*k4+(99/256).*k5+(1).*k7);
        
        sol5(:, j+1) = sol5(:, j) ...
                           +(31/384).*k1+(1125/2816).*k3+(9/32).*k4+(125/768).*k5+(5/66).*k6;
        sol6(:, j+1) = sol5(:, j) ...
                           +(7/1408).*k1+(1125/2816).*k3+(9/32).*k4+(125/768).*k5+(5/66).*k7+(5/66).*k8;
 
        %Time and stepsize updates
        time(j+1) = time(j) + stepsize;
        s = (tol*stepsize/(2*norm(sol6(:, j+1)-sol5(:, j+1)))).^(1/4);
        stepsize = min([s*stepsize max_stepsize]);
        j=j+1;
    end
        
end