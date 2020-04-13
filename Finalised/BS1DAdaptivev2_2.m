function [time, solution] = BS1DAdaptivev2_2(func, tspan, initial, H, kmax, atol, rtol)
    solution(:, 1) = initial;
    n_eqs = size(initial, 1);
    time(1) = tspan(1);
    T = zeros(kmax, kmax, n_eqs);
    S1=0.95;
    S2=0.95;
    kfinal=kmax;
    k = 1; 
    m = 1;
    % The issue is kfinal is undefined if it does all the step
    while time(m)<=tspan(2)
        err = 2; % Arbitrary large number bigger than 1
        while k <= kmax
            if err > 1 && k~=kmax
                for j=1:k
                    if j==1
                        % The first column of results
                        T(k, 1, :) = midpointMethod(func, time(m), solution(:,m), H, 2*k);
                    else
                        T(k, j, :) = T(k, j-1, :) + (T(k, j-1, :) - T(k-1, j-1, :))/((2*k/(2*(k-(j-1))))^2-1);
                    end 
                end
                if k ~= 1
                % See page 912 of numerical recipes of C, third edition
                    triangle = reshape(abs(T(k, k, :) - T(k, k-1, :)), [1, n_eqs]);
                    triangleDivScale = triangle/(atol+rtol*norm(reshape(T(k, k, :), [1, n_eqs])));
                    err = (n_eqs)^(-1/2)*norm(triangleDivScale);
                end
                k = k + 1;
            elseif err>1 && k==kmax
                H = H*S1*(S2/err)^(1/(2*k+1));
                k = 1;
            elseif err<=1 && err>0
                prev_H = H;
                H = H*S1*(S2/err)^(1/(2*k+1));
                kfinal = k-1;
                break
            elseif err==0
                prev_H = H;
                kfinal = k-1;
                break
            end
        end
        solution(:, m+1) = T(kfinal, kfinal, :);
        time(m+1) = time(m) + prev_H;
        m = m + 1;
    end
end