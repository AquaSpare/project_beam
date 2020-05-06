
close all

%exakt l�sning och begynnelsev�rde f�r fria r�nder

%konvergens f�r SBP4 

iter = 6;
N = 15;
errorProj = zeros(1,iter);
stepsProj = zeros(1,iter);
errorSAT = zeros(1,iter);
stepsSAT = zeros(1,iter);

x0 = 0;
xN = 1;
b = 1;
BC = 2;
T = 0.14;
order = 4;

exact = homo_beam_ana_2(BC);
for i = 1:iter
    [projSolution,tProj,xProj,hProj,kProj] = beam_homogenous((N*i)+1,x0,xN,T,b,order,BC,500,0.00001,0,2);
    stepsProj(i) = hProj;
    errorProj(i) = sqrt(hProj)*norm(exact(xProj,T)'-projSolution(:,end),2);
    [SATsolution, tSAT, xSAT, hSAT, kSAT] = SAT_enBalk((N*i)+1,x0,xN,T, b, order, BC, 500, 0.00001, 0, 2);
    stepsSAT(i) = hSAT;
    errorSAT(i) = sqrt(hSAT)*norm(exact(xSAT,T)'-SATsolution(:,end),2);
end

figure(1)
title('Homogenous beam, free boundary')
loglog(stepsProj,errorProj,'r*')
hold on
loglog(stepsSAT,errorSAT,'g*')
hold on
loglog(stepsProj,stepsProj.^4,'--')
hold off
legend('Projection', 'SAT','4th order slope')
xlabel('h')
ylabel('L2 norm of error')

%unders�ka felet �ver tid
                                
% [sol,t,x,h,k] = beam_homogenous(41,0,1,10,1,4,1,500,0.0001,0,1);
% errortime = zeros(1,length(t));
% 
% for i = 1:length(t)
%     errortime(1,i) = sqrt(h)*norm(exact(x,t(i)) - sol(:,i),2);        
% end
% 
% figure(2)
% 
% plot(t,errortime)
% xlabel('time')
% ylabel('L2 norm of error')
