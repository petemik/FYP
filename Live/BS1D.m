function next_step = BS1D(func, tspan, x0, H, kmax)
    initial = x0';
    n_eqs = size(initial, 1);
    total_steps = ceil((tspan(2)-tspan(1))/H);
    time = zeros(1, total_steps);
    time(1) = tspan(1);
    
    T = zeros(kmax, kmax, n_eqs);
    for k=1:kmax
        for j=1:k
            if j==1
                T(k, 1, :) = midpointMethod(func, time(1), initial, H, 2*k);
            else
                T(k, j, :) = T(k, j-1, :) + (T(k, j-1, :) - T(k-1, j-1, :))/((2*k/(2*(k-(j-1))))^2-1);
            end 
        end
    end
    next_step = T(kmax, kmax, :);
end