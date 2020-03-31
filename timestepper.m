function [solution, k] = timestepper(t0, T, h, A, u0, u0_t)

k = 0.1*h;

t = t0:k:T;
solution = zeros(length(u0),length(t));
solution(:,1) = u0;
solution(:,2) = k*u0_t' + 0.5*k^2*A*u0' + u0';
for i=3:length(t)
    solution(:,i) = k^2.*A*solution(:,i-1) + 2*solution(:,i-1) - solution(:,i-2);
end

