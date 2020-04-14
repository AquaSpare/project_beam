close all

%exakt l�sning och begynnelsev�rde f�r fria r�nder
u0 = @(x) cosh(1.50562*pi.*x) + cos(1.50562*pi.*x) - ((cos(1.50562*pi) -cosh(1.50562*pi))/(sin(1.50562*pi) - sinh(1.50562*pi)))*(sin(1.50562*pi.*x) +sinh(1.50562*pi.*x));
u_exact = @(x,t) real(exp(-1i*(22.3733)*t)*u0(x));

%konvergens f�r SBP4 

iter = 5;
N = 20;
errorProj = zeros(1,iter);
errorSAT = zeros(1,iter);
stepsProj = zeros(1,iter);
stepsSAT = zeros(1,iter);


for i = 1:iter
    [solProj,tProj,xProj,hProj,kProj] = beam_homogenous((N*i)+1,0,1,0.14,1,4,4,500,0.00001,0);
    exact = u_exact(xProj,tProj(end));
    stepsProj(i) = hProj;
    errorProj(i) = sqrt(1/(N*i))*norm(exact'-solProj(:,end),2);
    
    [solSAT, tSAT, xSAT, hSAT, kSAT] = SAT_enBalk((N*i)+1,0,1,0.14, 1, 4, 2, 500, 0.00001,0);
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