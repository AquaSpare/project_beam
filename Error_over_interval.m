%Error over interval to see where the errors are the highest.
N = 161;
x0 = -1;
xl=0;
xN = 1;

T = 0.14;
a1 = 1;
a2 = 100;
order = 4;
BC = 2;
k_ratio = 0.01;
timestepper = 2;

[u_exact_l, u_exact_r] = beam_ana(BC,a1,a2);

[w,t,x,h,k] = proj_projSAT_proj(N, x0, xl, xN, T, a1, a2, order, BC, 100, k_ratio, 1, 0, timestepper);
error_interval_projSAT = zeros(length(x),1);
u_exact = [u_exact_l(x0:h:xl,t(end)) u_exact_r(xl:h:xN,t(end))];
for i = 1:length(x)
    error_interval_projSAT(i) = norm(u_exact(i)-w(i,end),2);
end

% [w,t,x,h,k] = proj_projSAT2_proj(N, x0, xl, xN, T, a1, a2, order, BC, 100, k_ratio, 1, 0, timestepper);
% error_interval_projSAT2 = zeros(length(x),1);
% u_exact = [u_exact_l(x0:h:xl,t(end)) u_exact_r(xl:h:xN,t(end))];
% for i = 1:length(x)
%     error_interval_projSAT2(i) = norm(u_exact(i)-w(i,end),2);
% end
% 
% [w,t,x,h,k] = proj_projSAT3_proj(N, x0, xl, xN, T, a1, a2, order, BC, 100, k_ratio, 1, 0, timestepper);
% error_interval_projSAT3 = zeros(length(x),1);
% u_exact = [u_exact_l(x0:h:xl,t(end)) u_exact_r(xl:h:xN,t(end))];
% for i = 1:length(x)
%     error_interval_projSAT3(i) = norm(u_exact(i)-w(i,end),2);
% end

figure(1);
% plot(x,error_interval_projSAT,x,error_interval_projSAT2,x,error_interval_projSAT3);
plot(x,error_interval_projSAT);
title('Error over domain','FontSize',18);
xlabel('x','FontSize',16);
ylabel('L2 norm of error','FontSize',16);
hold on;
plot([0 0], [0.2e-6 0], 'r--', 'LineWidth',0.5);
% legend('ProjSAT','ProjSAT2','ProjSAT3','Inner boundary');
legend('ProjSAT','Inner boundary');
hold off;

figure(2);
plot(x0:h:xl,u_exact_l(x0:h:xl,length(t)*k-k),'r' , xl:h:xN, u_exact_r(xl:h:xN,length(t)*k-k),'r');
hold on;
plot(x, w(:,end),'*b');
title('solution at t=T');
