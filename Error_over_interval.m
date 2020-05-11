%Error over interval to see where the errors are the highest.
N = 41;
x0 = -1;
xl=0;
xN = 1;

T = 1000;
a1 = 1;
a2 = 4;
order = 4;
BC = 2;
k_ratio = 0.05;
timestepper = 2;

[u_exact_l, u_exact_r] = beam_ana(BC,a1,a2);
[w,t,x,h,k] = proj_projSAT_proj(N, x0, xl, xN, T, a1, a2, order, BC, 100, k_ratio, 1, 0, timestepper);

error_interval = zeros(length(x),1);
u_exact = [u_exact_l(x0:h:xl,length(t)*k-k) u_exact_r(xl:h:xN,length(t)*k-k)];
for i = 1:length(x)
    error_interval(i) = norm(u_exact(end,i)-w(i,end),2);
end
figure(1);
plot(x,error_interval)
title('error');

figure(2);
plot(x0:h:xl,u_exact_l(x0:h:xl,length(t)*k-k),'r' , xl:h:xN, u_exact_r(xl:h:xN,length(t)*k-k),'r');
hold on;
plot(x, w(:,end),'*b');
title('solution at t=T');
