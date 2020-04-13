function next_step = midpointMethod(f, tspan, initial, H, n)
%tspan really should be t0
% Stepsize is equivalent to big H in the numerical methods book
% n=1 for euler method
    h = H/n;
    n_eqs = size(initial, 1);
    z = zeros(n_eqs, n+1);
    t0 = tspan(1);
    z(:, 1) = initial;
    z(:, 2) = z(:, 1) + h*f(t0, z(:, 1));
    for m=2:n
         z(:, m+1) = z(:, m-1) + 2*h*f(t0+(m-1)*h, z(:, m));
    end
    next_step = 1/2*(z(:, n+1)+z(:, n)+h*f(t0+H, z(:, n+1)));
%The above is equivalent to y(x+H) in the book
end