function [solution, k, t, error] = timestepperv3(T, h, A, u0, u0_t, k_ratio,BC,a1,a2,x0,xl,xN)

k = k_ratio*h^2;

t = 0:k:T;
solution = zeros(length(u0),3);
solution(:,1) = u0;
solution(:,2) = k*u0_t' + 0.5*k^2*A*u0' + u0';
[u_exact_l, u_exact_r] = beam_ana(BC,a1,a2);

saves = 2000;
error = zeros(1,ceil(length(t)/saves));

j = 0;
for i=3:length(t)
    solution(:,3) = k^2.*A*solution(:,2) + 2*solution(:,2) - solution(:,1);
    solution(:,1) = solution(:,2);
    solution(:,2) = solution(:,3);
    if mod(i,saves) == 0
        j = j+1;
        error(j) = sqrt(h)*norm([u_exact_l(x0:h:xl,t(i)) u_exact_r(xl:h:xN,t(i))]-solution(:,3)',2);
    end
end