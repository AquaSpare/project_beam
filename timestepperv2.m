function [solution, k, t] = timestepperv2(T, h, A, u0, u0_t, k_ratio,order)

%k = (h^2)/sqrt(max(abs(eig(h^4.*A))));
k = k_ratio*h^2;

t = 0:k:T;
solution = zeros(length(u0),3);

if order == 2
    solution(:,1) = u0;
    solution(:,2) = k*u0_t' + 0.5*k^2*A*u0' + u0';
    for i=3:length(t)
        solution(:,3) = k^2.*A*solution(:,2) + 2*solution(:,2) - solution(:,1);
        solution(:,1) = solution(:,2);
        solution(:,2) = solution(:,3);
    end
elseif order == 4
    solution(:,1) = u0;
    solution(:,2) =  k*u0_t' + (k^3/6).*A*u0_t' + 0.5*k^2*A*u0' + u0';
    for i=3:length(t)
        solution(:,3) = (k^4/12)*A^2*solution(:,2)+k^2.*A*solution(:,2) + 2*solution(:,2) - solution(:,1);
        solution(:,1) = solution(:,2);
        solution(:,2) = solution(:,3);
    end
end