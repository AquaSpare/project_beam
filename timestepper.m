function [solution, k, t] = timestepper(T, h, A, u0, u0_t, k_ratio)

%k = 2/sqrt(max(abs(eig(h^2.*A))));
k = k_ratio*h;

t = 0:k:T;
solution = zeros(length(u0),length(t));
solution(:,1) = u0;
solution(:,2) = k*u0_t' + 0.5*k^2*A*u0' + u0';
for i=3:length(t)
    solution(:,i) = k^2.*A*solution(:,i-1) + 2*solution(:,i-1) - solution(:,i-2);
end

