function next_step = BS1DAdaptive(func, tspan, x0, H, kmax, atol, rtol)
    initial = x0';
    n_eqs = size(initial, 1);
    total_steps = ceil((tspan(2)-tspan(1))/H);
    time = zeros(1, total_steps);
    time(1) = tspan(1);
    err = 2; % Arbitrary large number bigger than 1
    T = zeros(kmax, kmax, n_eqs);
    % The issue is kfinal is undefined if it does all the step 
    for k=1:kmax
        if err > 1
            for j=1:k
                if j==1
                    % The first column of results
                    T(k, 1, :) = midpointMethod(func, time(1), initial, H, 2*k);
                else
                    T(k, j, :) = T(k, j-1, :) + (T(k, j-1, :) - T(k-1, j-1, :))/((2*k/(2*(k-(j-1))))^2-1);
                end 
            end
        else 
           kfinal = k-1;
           break
        end
        if k~= 1
        % See page 912 of numerical recipes of C, third edition
            triangle = reshape(abs(T(k, k, :) - T(k, k-1, :)), [1, 36]);
            triangleDivScale = triangle/(atol+rtol*norm(reshape(T(k, k, :), [1, 36])));
            err = (n_eqs)^(-1/2)*norm(triangleDivScale);
        end
    end   
    
    next_step = T(kfinal, kfinal, :);
end