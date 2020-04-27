%%% J�mf�r de 3 olika metoderna. �ndra N f�r antal punkter och iter f�r
%%% sista antalet punkter ( iter*N ) 

iter = 7;
N = 15;
errorProj = zeros(1,iter);
errorSAT = zeros(1,iter);
stepsProj = zeros(1,iter);
stepsSAT = zeros(1,iter);
stepsProjSAT = zeros(1,iter);
errorProjSAT = zeros(1,iter);

BC = 1;
a1 = 1;
a2 = 4;
T = 0.14;
order = 4;

x0 = -1;
xl = 0;
xN = 1;

[u,v] = beam_ana(BC,a1,a2);

for i = 1:iter
    [solProj,tProj,xProj,hProj,kProj] = proj_proj_proj((N*i)+1, x0, xl, xN, T, a1, a2, order, BC, 10, 0.0001, 1, 0, 2);
    exact = [u(x0:hProj:xl,tProj(end)) v(xl:hProj:xN,tProj(end))];
    stepsProj(i) = hProj;
    errorProj(i) = sqrt(hProj)*norm(exact'-solProj(:,end),2);
                                        
%     [solSAT, tSAT, xSAT, hSAT, kSAT] = SAT_projSAT_SAT((N*i)+1,x0,xl,xN,T, a1, a2, order, BC, 10,0.0001,1,0,2);
%     exact = [u(x0:hSAT:xl,tSAT(end)) v(xl:hSAT:xN,tSAT(end))];
%     stepsSAT(i) = hSAT;
%     errorSAT(i) = sqrt(hSAT)*norm(exact'-solSAT(:,end),2);
    
    [solProjSAT, tProjSAT, xProjSAT, hProjSAT, kProjSAT] = proj_projSAT_proj((N*i)+1,x0,xl,xN,T, a1, a2, order, BC, 100,0.0001,1,0,2);
    exact = [u(x0:hProjSAT:xl,tProjSAT(end)) v(xl:hProjSAT:xN,tProjSAT(end))];
    stepsProjSAT(i) = hProjSAT;
    errorProjSAT(i) = sqrt(hProjSAT)*norm(exact'-solProjSAT(:,end),2);
end
figure(1)
title('Inner boundary, clamped ends')
loglog(stepsProj,errorProj,'r*')
hold on
% loglog(stepsSAT,errorSAT,'b*')
hold on
loglog(stepsProjSAT,errorProjSAT, 'g*')
xlabel('log(h)')
ylabel('L2 norm of error')
legend('Pure Projection', 'Outer Projection, inner Mix');