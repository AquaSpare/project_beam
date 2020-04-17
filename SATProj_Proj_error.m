%%% Jämför de 3 olika metoderna. Ändra N för antal punkter och iter för
%%% sista antalet punkter ( iter*N ) 

iter = 7;
N = 15;
errorProj = zeros(1,iter);
errorSAT = zeros(1,iter);
stepsProj = zeros(1,iter);
stepsSAT = zeros(1,iter);
stepsProjSAT = zeros(1,iter);
errorProjSAT = zeros(1,iter);

[u,v] = beam_ana();

for i = 1:iter
    [solProj,tProj,xProj,hProj,kProj] = proj_proj_proj((N*i)+1,-1,0,1,0.14,1,4,4,2,100,0.0001,1,0,2);
    exact = [u(xProj(1:(N*i)+1),tProj(end)) v(xProj((N*i)+2:end),tProj(end))];
    stepsProj(i) = hProj;
    errorProj(i) = sqrt(1/(N*i))*norm(exact'-solProj(:,end),2);
                                        
    [solSAT, tSAT, xSAT, hSAT, kSAT] = SAT_projSAT_SAT((N*i)+1,-1,0,1,0.14, 1, 4, 4, 2, 100,0.0001,1,0,2);
    exact = [u(xSAT(1:(N*i)+1),tSAT(end)) v(xSAT((N*i)+2:end),tSAT(end))];
    stepsSAT(i) = hSAT;
    errorSAT(i) = sqrt(1/(N*i))*norm(exact'-solSAT(:,end),2);
    
    [solProjSAT, tProjSAT, xProjSAT, hProjSAT, kProjSAT] = proj_projSAT_proj((N*i)+1,-1,0,1,0.14, 1, 4, 4, 2, 100,0.0001,1,0,2);
    exact = [u(xProjSAT(1:(N*i)+1),tProjSAT(end)) v(xProjSAT((N*i)+2:end),tProjSAT(end))];
    stepsProjSAT(i) = hProjSAT;
    errorProjSAT(i) = sqrt(1/(N*i))*norm(exact'-solProjSAT(:,end),2);
end
figure(1)
loglog(stepsProj,errorProj,'r*')
hold on
loglog(stepsSAT,errorSAT,'b*')
hold on
loglog(stepsProjSAT,errorProjSAT, 'g*')
xlabel('log(h)')
ylabel('L2 norm of error')
legend('Pure Projection', 'Outer SAT, inner Mix', 'Outer Projection, inner Mix');