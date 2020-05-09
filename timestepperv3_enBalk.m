function [solution, k, t, error] = timestepperv3_enBalk(T, h, A, u0, u0_t, k_ratio,BC,a,x0,xN)

k = k_ratio*h^2;

t = 0:k:T;
solution = zeros(length(u0),3);
solution(:,1) = u0;
solution(:,2) = k*u0_t' + 0.5*k^2*A*u0' + u0';
u_exact = homo_beam_ana_2(a,BC);

saves = 2000;
error = zeros(1,ceil(length(t)/saves));

j = 0;
for i=3:length(t)
    solution(:,3) = k^2.*A*solution(:,2) + 2*solution(:,2) - solution(:,1);
    solution(:,1) = solution(:,2);
    solution(:,2) = solution(:,3);
    if mod(i,saves) == 0
        j = j+1;
        error(j) = sqrt(h)*norm(u_exact(x0:h:xN,t(i))-solution(:,3)',2);
    end
end