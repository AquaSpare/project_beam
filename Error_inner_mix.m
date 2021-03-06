%%% Jämför de 3 olika metoderna. Ändra N för antal punkter och iter för
%%% sista antalet punkter ( iter*N ) 
close all
iter = 5;
N = 30;
errorProjSAT = zeros(1,iter);
errorProjSAT2 = zeros(1,iter);
errorProjSAT3 = zeros(1,iter);
errorProj = zeros(1,iter);
stepsProj = zeros(1,iter);
steps = zeros(1,iter);

BC = 1;
a1 = 1;
a2 = 100;
T = 0.14;
order = 2;
kratio = 0.01;
timestepperOrder = 4;

x0 = -1;
xl = 0;
xN = 1;

[u,v] = beam_ana(BC,a1,a2);

for i = 1:iter
    [sol,t,x,h,k] = proj_projSAT_proj((N*i)+1, x0, xl, xN, T, a1, a2, order, BC, 100, kratio, 1, 0, 2,timestepperOrder);
    exact = [u(x0:h:xl,t(end)) v(xl:h:xN,t(end))];
    steps(i) = h;
    errorProjSAT(i) = sqrt(h)*norm(exact'-sol(:,end),2);
    
    [sol, t, x, h, k] = proj_projSAT2_proj((N*i)+1, x0, xl, xN, T, a1, a2, order, BC, 100, kratio, 1, 0, 2,timestepperOrder);
    
    errorProjSAT2(i) = sqrt(h)*norm(exact'-sol(:,end),2);
    
    [sol, t, x, h, k] = proj_projSAT3_proj((N*i)+1, x0, xl, xN, T, a1, a2, order, BC, 100, kratio, 1, 0, 2,timestepperOrder);
    
    errorProjSAT3(i) = sqrt(h)*norm(exact'-sol(:,end),2);
    
    [sol, t, x, h, k] = proj_proj_proj((N*i)+1, x0, xl, xN, T, a1, a2, order, BC, 100, kratio, 1, 0, 2,timestepperOrder);
    
    errorProj(i) = sqrt(h)*norm(exact'-sol(:,end),2);
    
end

figure(1)
loglog(steps,errorProj, 'm^')
hold on
loglog(steps,errorProjSAT2, 'rx')
loglog(steps,errorProjSAT3, 'go')
loglog(steps,errorProjSAT,'b*')
loglog(steps,100*steps.^2,'k--')
hold off
title('Inner boundary, clamped outer','FontSize',18);
legend('Pure projection','2nd derivative SAT','3rd derivative SAT','2nd & 3rd derivative SAT','2th order convergence','Location', 'northwest');
xlabel('h','FontSize',16);
ylabel('L2 norm of error','FontSize',16)
