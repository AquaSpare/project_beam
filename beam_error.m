
close all


iter = 8;
N = 10;
error = zeros(1,iter);
steps = zeros(1,iter);
BC = 2;
a1 = 1;
a2 = 4;
T = 0.5;
order = 4;

x0 = -1;
xl = 0;
xN = 1;

[u_exact_l, u_exact_r] = beam_ana(BC,1,4);

for i = 1:iter
    [solution,t,x,h,k] = proj_proj_proj((N*i)+1, x0, xl, xN, T, a1, a2, order, BC, 10, 0.0001, 1, 0, 2);
    exact = [u_exact_l(x0:h:xl,t(end)) u_exact_r((xl:h:xN),t(end))];
    steps(i) = h;
    error(i) = sqrt(h)*norm(exact'-solution(:,end),2);
end
figure(1)
loglog(steps,error,'r*')
xlabel('log(h)')
ylabel('L2 norm of error')



%tar fram felet över tid, ser hur det växer.
% T = 1;
% [solution,t,x,h,k] = proj_proj_proj(41, x0, xl, xN, T, a1, a2, order, BC, 10, 0.0001, 1, 0, 1);
% exact = zeros(size(solution));
% errortime = zeros(1,length(t));
% for i =1:length(t)
%     exact = [u_exact_l(x0:h:xl,t(i)) u_exact_r(xl:h:xN,t(i))];
%     errortime(i) = sqrt(h)*norm(exact' - solution(:,i),2);
% end
% figure(2)
% plot(t,errortime)
% xlabel('time')
% ylabel('L2 norm of error')

