function [solution, k, t] = timestepper(t0, T, h, A, u0, u0_t)

%k = 2*h/sqrt(max(abs(eig(A))));
<<<<<<< HEAD
 k = 0.0001*h;
=======
 k = 0.001*h;
>>>>>>> b37fb26f351043c6734435c6c131adcb39491323

t = t0:k:T;
solution = zeros(length(u0),length(t));
solution(:,1) = u0;
solution(:,2) = k*u0_t' + 0.5*k^2*A*u0' + u0';
for i=3:length(t)
    solution(:,i) = k^2.*A*solution(:,i-1) + 2*solution(:,i-1) - solution(:,i-2);
end

