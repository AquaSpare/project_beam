%Wavefunction u_tt = c^2u_xx
close all
pause on 
%%%Domain%%%
N = 161;
x0 = -1;
xl = 0;
xN = 1;
h = (xN-xl)/(N-1);
x1 = x0:h:xl;
x2 = xl:h:xN;
x = [x1 x2];
t0 = 0;
T = 10;

%%%%%%%%%%%%%

%%%Constant%%%
c1 = 1;
c2 = 2;
% lambda = pi*c1/x2;
%%%%%%%%%%%

%%%Boundary conditions%%%
% u_x0 = 0; u_xN = 0; % Neumann

%%%Initial data%%%
% u0 = @(x) sin(lambda.*x/c) + 2*sin(2*lambda.*x/c) + sin(3*lambda.*x/c) + sin(4*lambda.*x/c) + sin(5*lambda.*x/c) + sin(6*lambda.*x/c);
x0 = -1/4;
r0 = 1/30;
u0 = @(x) 5*exp(-(x0-x).^2/r0^2);
u0_t = 0;

%%%Exact solution%%%
u_exact = @(x,t) cos(lambda.*t)*sin(lambda.*x/c) + cos(2*lambda.*t)*2*sin(2*lambda.*x/c) + cos(3*lambda.*t)*sin(3*lambda.*x/c)+ cos(4*lambda.*t)*sin(4*lambda.*x/c) + cos(5*lambda.*t)*sin(5*lambda.*x/c) + cos(6*lambda.*t)*sin(6*lambda.*x/c);
[D2, H, HI, M, e1, eN, d1, dN] = SBP4(N, h);

L = [eN -e1; c1^2*dN -c2^2*d1; e1 zeros(1,N); zeros(1,N) eN];
H = [H zeros(N); zeros(N) H];
HI = inv(H);

P = [eye(N) zeros(N); zeros(N) eye(N)] - HI*L'*inv(L*HI*L')*L;

A = (P*[c1^2*D2 zeros(N); zeros(N) c2^2*D2]*P);

w0 = [u0(x1) u0(x2)];
w0_t = 0;

[w,k,t] = timestepper(t0,T,h,A,w0,w0_t);
figure(1);
for i = 1:1:length(t)
    plot(x,w(:,i), '*b');
%     hold on;
%     plot(x,u_exact(x,t(i)), 'r');
%     hold off;
    axis([-1 1 -6 6]);
    pause(0.00000001);
end
