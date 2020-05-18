%%% Jämför de 3 olika metoderna. Ändra N för antal punkter och iter för
%%% sista antalet punkter ( iter*N ) 
close all
iter = 5;
N = 20;
error2 = zeros(1,iter);
error4 = zeros(1,iter);
error6 = zeros(1,iter);
time2 = zeros(1,iter);
time4 = zeros(1,iter);
time6 = zeros(1,iter);

BC = 1;
a1 = 1;
a2 = 4;
T = 3;
timestepperOrder = 4;

x0 = -1;
xl = 0;
xN = 1;

[u,v] = beam_ana(BC,a1,a2);

for i = 1:iter
    tic;
    [sol,t,x,h,k] = proj_projSAT_proj((N*i)+1, x0, xl, xN, T, a1, a2, 2, BC, 100, 0, 1, 0, 2,timestepperOrder);
    time2(i) = toc;
    exact = [u(x0:h:xl,t(end)) v(xl:h:xN,t(end))];
    error2(i) = sqrt(h)*norm(exact'-sol(:,end),2);
    
    tic;
    [sol,t,x,h,k] = proj_projSAT_proj((N*i)+1, x0, xl, xN, T, a1, a2, 4, BC, 100, 0, 1, 0, 2,timestepperOrder);
    time4(i) = toc;
    exact = [u(x0:h:xl,t(end)) v(xl:h:xN,t(end))];
    error4(i) = sqrt(h)*norm(exact'-sol(:,end),2);
    
    tic;
    [sol,t,x,h,k] = proj_projSAT_proj((N*i)+1, x0, xl, xN, T, a1, a2, 6, BC, 100, 0, 1, 0, 2,timestepperOrder);
    time6(i) = toc;
    exact = [u(x0:h:xl,t(end)) v(xl:h:xN,t(end))];
    error6(i) = sqrt(h)*norm(exact'-sol(:,end),2);
    
end

figure(1)
loglog(time2,error2, 'ro')
hold on
loglog(time4,error4, 'bs')
loglog(time6,error6, 'g^')
hold off
title('Runtime, inner boundary','FontSize',18);
legend('Order 2','Order 4','Order 6','Location', 'northeast');
xlabel('runtime','FontSize',16);
ylabel('L2 norm of error','FontSize',16)
