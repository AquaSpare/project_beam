function [solution, k, t] = timestepperv2(t0, T, h, A, u0, u0_t)

k = 2*h/sqrt(max(abs(eig(A))));
%k = 0.1*h;

t = t0:k:T;
solution = zeros(length(u0),3);
solution(:,1) = u0;
solution(:,2) = k*u0_t' + 0.5*k^2*A*u0' + u0';
for i=3:length(t)
    solution(:,3) = k^2.*A*solution(:,2) + 2*solution(:,2) - solution(:,1);
    solution(:,1) = solution(:,2);
    solution(:,2) = solution(:,3);
end