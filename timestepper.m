function [u, k, t] = timestepper(T, h, A, u0, u0_t, k_ratio)

% k = 2*h/sqrt(max(abs(eig(A))));
k = k_ratio*h;

t = 0:k:T;
u = zeros(length(u0),length(t));
u(:,1) = u0;
u(:,2) = k*u0_t' + 0.5*k^2*A*u0' + u0';
for i=3:length(t)
    u(:,i) = k^2.*A*u(:,i-1) + 2*u(:,i-1) - u(:,i-2);
end

