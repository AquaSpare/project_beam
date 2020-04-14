
close all

%exakt lösning och begynnelsevärde för fria ränder
u0 = @(x) cosh(1.50562*pi.*x) + cos(1.50562*pi.*x) - ((cos(1.50562*pi) -cosh(1.50562*pi))/(sin(1.50562*pi) - sinh(1.50562*pi)))*(sin(1.50562*pi.*x) +sinh(1.50562*pi.*x));
u_exact = @(x,t) real(exp(-1i*(22.3733)*t)*u0(x));

%konvergens för SBP4 

iter = 5;
N = 10;
error = zeros(1,iter);
steps = zeros(1,iter);


for i = 1:iter
    [sol,t,x,h,k] = beam_homogenous((N*i)+1,0,1,0.2,1,4,4,500,0.00001,0);
    exact = u_exact(x,t(end));
    steps(i) = h;
    error(i) = sqrt(1/(N*i))*norm(exact-sol(:,end),2);
end

figure(1)
loglog(steps,error,'r*')
xlabel('log(h)')
ylabel('L2 norm of error')

%undersöka felet över tid

[sol,t,x,h,k] = beam_homogenous(41,0,1,10,1,6,4,500,0.0001,0);
exact = zeros(size(sol));
errortime = zeros(1,length(t));

for i = 1:length(t)
    exact(:,i) = u_exact(x,t(i));
    errortime(1,i) = sqrt(1/40)*norm(exact(:,i) - sol(:,i),2);        
end

figure(2)

plot(t,errortime)
xlabel('time')
ylabel('L2 norm of error')
