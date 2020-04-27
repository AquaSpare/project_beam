%Error over interval to see where the errors are the highest.
N = 81;
x0 = -1;
xl=0;
xN = 1;

[u_exact_l, u_exact_r] = beam_ana(2,1,4);
[w,t,x,h,k] = proj_projSAT_proj(N, x0, xl, xN, 0.1, 1, 4, 4, 2, 100, 0.0001, 1, 0, 1);

error_interval = zeros(length(x),1);
u_exact = [u_exact_l(x0:h:xl,0.1) u_exact_r(xl:h:xN,0.1)];
for i = 1:length(x)
    error_interval(i) = norm(u_exact(end,i)-w(i,end),2);
end
plot(x,error_interval)
