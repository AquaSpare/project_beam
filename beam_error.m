%[w,t,x,h,k] = beam_eq_projection(N, x0, xl, xN, T, b1, b2, order, BC, plotspeed, k_ratio)

close all

%Exakt lösning (typ) för DBE med fria ränder och längd 1.
u0 = @(x) cosh(1.50562*pi.*x) + cos(1.50562*pi.*x) - ((cos(1.50562*pi) -cosh(1.50562*pi))/(sin(1.50562*pi) - sinh(1.50562*pi)))*(sin(1.50562*pi.*x) +sinh(1.50562*pi.*x));
u_exact = @(x,t) real(exp(-1i*(22.3733)*t)*u0(x));

%Konvergensstudie för SBP-order 4
iter = 8;
N = 10;
error = zeros(1,iter);
steps = zeros(1,iter);
for i = 1:iter
    [solution,t,x,h,k] = beam_eq_projection((N*i)+1,0,0.5,1,0.18,1,1,4,2,2000,0.00001);
    exact = u_exact(x,t(end))';
    steps(i) = h;
    error(i) = sqrt(1/(N*i))*norm(exact-solution(:,end),2);
end
figure(1)
loglog(steps,error,'r*')
xlabel('log(h)')
ylabel('L2 norm of error')



% tar fram felet över tid, ser hur det växer.
[solution,t,x,h,k] = beam_eq_projection_dav(41,0,0.5,1,0,20,1,1,4,2);
exact = zeros(size(solution));
errortime = zeros(1,length(t));
for i =1:length(t)
    exact(:,i) = u_exact(x,t(i));
    errortime(i) = sqrt(1/40)*norm(exact(:,i) - solution(:,i),2);
end
figure(2)
plot(t,errortime)
xlabel('time')
ylabel('L2 norm of error')

