%function [w,t,x,h,k] = beam_eq_projection_dav(N,x0,xl,xN,t0,T,b1,b2,order,BC)

close all
u0 = @(x) cosh(1.50562*pi.*x) + cos(1.50562*pi.*x) - ((cos(1.50562*pi) -cosh(1.50562*pi))/(sin(1.50562*pi) - sinh(1.50562*pi)))*(sin(1.50562*pi.*x) +sinh(1.50562*pi.*x));
u_exact = @(x,t) real(exp(-1i*(22.3733)*t)*u0(x));
% [solution51,t41,x41] = beam_eq_projection_dav(41,0,0.5,1,0,0.1,1,1,4,2);
% [solution81,t81,x81] = beam_eq_projection_dav(81,0,0.5,1,0,0.1,1,1,4,2);
% [solution161,t5161,x161] = beam_eq_projection_dav(161,0,0.5,1,0,0.1,1,1,4,2);
iter = 6;
N = 20;
error = zeros(1,iter);
steps = zeros(1,iter);
% for i = 1:iter
%     [solution,t,x,h,k] = beam_eq_projection_dav((N*i)+1,0,0.5,1,0,0.1,1,1,4,2);
%     exact = u_exact(x,t(end))';
%     steps(i) = h;
%     error(i) = sqrt(1/(N*i))*norm(exact-solution(:,end),2);
% end
% loglog(steps,error,'*')
% axis([0 0.03 10^-6 10^-2])

[solution,t,x,h,k] = beam_eq_projection_dav(41,0,0.5,1,0,2,1,1,4,2);
exact = zeros(size(solution));
errortime = zeros(1,length(t));
for i =1:length(t)
    exact(:,i) = u_exact(x,t(i));
    errortime(i) = sqrt(1/100)*norm(exact(:,i) - solution(:,i),2);
end


