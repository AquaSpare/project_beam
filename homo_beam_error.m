
close all

%exakt l�sning och begynnelsev�rde f�r fria r�nder
u0 = @(x) cosh(1.50562*pi.*x) + cos(1.50562*pi.*x) - ((cos(1.50562*pi) -cosh(1.50562*pi))/(sin(1.50562*pi) - sinh(1.50562*pi)))*(sin(1.50562*pi.*x) +sinh(1.50562*pi.*x));
u_exact = @(x,t) real(exp(-1i*(22.3733)*t)*u0(x));

%konvergens f�r SBP4 

iter = 6;
N = 10;
error = zeros(1,iter);
steps = zeros(1,iter);

figure(1)
for i = 1:iter
    [sol,t,x,h,k] = beam_homogenous((N*(i+1))+1,0,1,0.2,1,4,4,500,0.00001,0);
    exact = u_exact(x,t(end));
    steps(i) = h;
    error(i) = sqrt(1/(N*i))*norm(exact-sol(:,end),2);
end

figure(2)
loglog(steps,error,'r*')
xlabel('log(h)')
ylabel('L2 norm of error')