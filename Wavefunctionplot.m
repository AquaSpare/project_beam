%Wavefunction u_tt = c^2u_xx
close all
pause on 
%%%Domain%%%
N = 201;
x0 = 0;
xN = 1;
h = (xN-x0)/(N-1);
x = x0:h:xN;
t0 = 0;
T = 1;

%%%%%%%%%%%%%

%%%Constant%%%
c = 1;
lambda = pi*c/xN;
%%%%%%%%%%%

%%%Boundary conditions%%%
u_x0 = 0; u_xN = 0; % Neumann

%%%Initial data%%%
u0 = cos(lambda.*x./c);
u0_t = 0;

%%%Exact solution%%%
u_exact = @(x,t) cos(lambda.*t).*cos(lambda*x/c);

[D2, H, HI, e1, eN, d1, dN] = SBP4(N, h);


A = c^2.*(D2 + HI*(e1'*d1)- HI*(eN'*dN));

[u,k] = timestepper(t0,T,h,A,u0,u0_t);
figure(1);
for i = 1:length(u)
    plot(x,u(:,i));
    pause(0.0001);
end

