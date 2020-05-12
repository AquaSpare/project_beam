
close all

%exakt l�sning och begynnelsev�rde f�r fria r�nder

%konvergens f�r SBP4 

iter = 4;
N = 20;
errorProj = zeros(1,iter);
stepsProj = zeros(1,iter);
errorSAT = zeros(1,iter);
stepsSAT = zeros(1,iter);

x0 = 0;
xN = 1;
a = 1;
BC = 3;
T = 10;
order = 4;
kratio = 0.1;

exact = homo_beam_ana_2(a,BC);
for i = 1:iter
    [projSolution,tProj,xProj,hProj,kProj] = beam_homogenous((N*i)+1,x0,xN,T,a,order,BC,500, kratio,1,0,2);
    stepsProj(i) = hProj;
    errorProj(i) = sqrt(hProj)*norm(exact(xProj,tProj(end))'-projSolution(:,end),2);
    [SATsolution, tSAT, xSAT, hSAT, kSAT] = SAT_enBalk((N*i)+1,x0,xN,T, a, order, BC, 500, kratio, 1, 0, 2);
    stepsSAT(i) = hSAT;
    errorSAT(i) = sqrt(hSAT)*norm(exact(xSAT,tProj(end))'-SATsolution(:,end),2);
end

figure(1)

loglog(stepsProj,errorProj,'*')
hold on
loglog(stepsSAT,errorSAT,'x')
hold off
if BC == 1 
    title('Homogenous beam, clamped boundary','FontSize',18)
elseif BC == 2
    title('Homogenous beam, free boundary','FontSize',18)
end
legend('Projection', 'SAT','Location','northwest')
xlabel('h','FontSize',16)
ylabel('L2 norm of error','FontSize',16)

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
