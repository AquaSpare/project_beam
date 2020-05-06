%%% Jämför de 3 olika metoderna. Ändra N för antal punkter och iter för
%%% sista antalet punkter ( iter*N ) 
close all
iter = 5;
N = 40;
errorProj = zeros(1,iter);
errorSAT = zeros(1,iter);
stepsProj = zeros(1,iter);
stepsSAT = zeros(1,iter);
stepsProjSAT = zeros(1,iter);
errorProjSAT = zeros(1,iter);
stepsSAT_proj_SAT = zeros(1,iter);
errorSAT_proj_SAT = zeros(1,iter);

BC = 2;
a1 = 1;
a2 = 4;
T = 0.14;
order = 4;

x0 = -1;
xl = 0;
xN = 1;

[u,v] = beam_ana(BC,a1,a2);

for i = 1:iter
    if BC == 1
    [solProj,tProj,xProj,hProj,kProj] = proj_proj_proj((N*i)+1, x0, xl, xN, T, a1, a2, order, BC, 10, 0.0001, 1, 0, 2);
    exact = [u(x0:hProj:xl,tProj(end)) v(xl:hProj:xN,tProj(end))];
    stepsProj(i) = hProj;
    errorProj(i) = sqrt(hProj)*norm(exact'-solProj(:,end),2);
    
    [solProjSAT, tProjSAT, xProjSAT, hProjSAT, kProjSAT] = proj_projSAT_proj((N*i)+1,x0,xl,xN,T, a1, a2, order, BC, 100,0.0001,1,0,2);
    exact = [u(x0:hProjSAT:xl,tProjSAT(end)) v(xl:hProjSAT:xN,tProjSAT(end))];
    stepsProjSAT(i) = hProjSAT;
    errorProjSAT(i) = sqrt(hProjSAT)*norm(exact'-solProjSAT(:,end),2);
    
    elseif BC == 2                                 
    [solSAT_projSAT_SAT, tSAT_projSAT_SAT, xSAT_projSAT_SAT, hSAT_projSAT_SAT, kSAT_projSAT_SAT] = SAT_projSAT_SAT((N*i)+1,x0,xl,xN,T, a1, a2, order, BC, 10,0.0001,1,0,2);
    exact = [u(x0:hSAT_projSAT_SAT:xl,tSAT_projSAT_SAT(end)) v(xl:hSAT_projSAT_SAT:xN,tSAT_projSAT_SAT(end))];
    stepsSAT(i) = hSAT_projSAT_SAT;
    errorSAT(i) = sqrt(hSAT_projSAT_SAT)*norm(exact'-solSAT_projSAT_SAT(:,end),2);
    
    [solSAT_proj_SAT, tSAT_proj_SAT, xSAT_proj_SAT, hSAT_proj_SAT, kSAT_proj_SAT] = SAT_proj_SAT((N*i)+1,x0,xl,xN,T, a1, a2, order, BC, 100,0.0001,1,0,2);
    exact = [u(x0:hSAT_proj_SAT:xl,tSAT_proj_SAT(end)) v(xl:hSAT_proj_SAT:xN,tSAT_proj_SAT(end))];
    stepsSAT_proj_SAT(i) = hSAT_proj_SAT;
    errorSAT_proj_SAT(i) = sqrt(hSAT_proj_SAT)*norm(exact'-solSAT_proj_SAT(:,end),2);
    end

end
figure(1)
if BC == 1
    loglog(stepsProj,errorProj,'b*')
    hold on
    loglog(stepsProjSAT,errorProjSAT, 'rx')
    hold off
elseif BC == 2
    loglog(stepsSAT_proj_SAT, errorSAT_proj_SAT,'b*');
    hold on
    loglog(stepsSAT,errorSAT,'rx')
    hold off
end

if BC == 1
    title('Inner boundary, clamped ends','FontSize',18)
    legend('Pure Projection','Outer Projection, inner mix','Location','northwest')
elseif BC == 2
    title('Inner boundary, free ends','FontSize',18)
    legend('Outer SAT, inner Projection', 'Outer SAT, inner mix','Location','northwest')
end
xlabel('h','FontSize',16);
ylabel('L2 norm of error','FontSize',16)
