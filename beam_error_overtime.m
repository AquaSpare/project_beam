clear all
close all
pause on

T = 5;
b = 1;
N = 21;
order = 6;
BC = 1;
kratio = 0.0001;
x0 = 0;
xN = 1;

[u{1,1}, t, x, h, k] = SAT_enBalk(N,x0,xN,T,b,order, BC, 500, kratio,0, 1);
[u{1,2},t,x,h,k] = beam_homogenous(N,x0,xN,T,b,order,BC,500,kratio,0, 1);
[u_exact,u0,u0_t] = homo_beam_ana(x,BC);
u_exact_sol = zeros(size(t));

errortime = zeros(2,length(t));

for i = 1:length(t)
    u_exact_sol = u_exact(x,t(i));
    errortime(1,i) = sqrt(1/N)*norm(u_exact_sol' - u{1,1}(:,i),2);
    errortime(2,i) = sqrt(1/N)*norm(u_exact_sol' - u{1,2}(:,i),2);
    i;
end

% plot(t,errortime(1,:),'b')
% hold on
% plot(t,errortime(2,:),'r')
% 
% legend('SAT','projection')
% xlabel('Time (s)')
% ylabel('L2 norm of error')
