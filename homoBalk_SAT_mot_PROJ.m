close all
Beta = 4.73004074486270402602404810;
%exakt l�sning och begynnelsev�rde f�r fria r�nder
u0 = @(x) cosh(Beta.*x) + cos(Beta.*x) - ((cos(Beta) -cosh(Beta))/(sin(Beta) - sinh(Beta)))*(sin(Beta.*x) +sinh(Beta.*x));
u_exact = @(x,t) real(exp(-1i*(Beta^2)*t)*u0(x));

%konvergens f�r SBP4 

iter = 1;
N = 400;
errorProj = zeros(1,iter);
errorSAT = zeros(1,iter);
stepsProj = zeros(1,iter);
stepsSAT = zeros(1,iter);

for i = 1:iter
    [solProj,tProj,xProj,hProj,kProj] = beam_homogenous((N*i)+1,0,1,0.14,1,4,4,500,0.0001,0,2);
    exact = u_exact(xProj,tProj(end));
    stepsProj(i) = hProj;
    errorProj(i) = sqrt(1/(N*i))*norm(exact'-solProj(:,end),2);
    
    [solSAT, tSAT, xSAT, hSAT, kSAT] = SAT_enBalk((N*i)+1,0,1,0.14, 1, 4, 2, 500, 0.0001,0,2);
    exact = u_exact(xSAT,tSAT(end));
    stepsSAT(i) = hSAT;
    errorSAT(i) = sqrt(1/(N*i))*norm(exact'-solSAT(:,end),2);
end

figure(1)
loglog(stepsProj,errorProj,'r*')
hold on
loglog(stepsSAT,errorSAT,'b*')
xlabel('log(h)')
ylabel('L2 norm of error')
legend('Projection', 'SAT');