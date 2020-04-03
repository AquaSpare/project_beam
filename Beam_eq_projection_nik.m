%Beamequation u_tt = -bu_xxxx
close all
pause on

%%%Domain%%%
N = 81;
x0 = -1;
xl = 0;
xN = 1;
h = (xN-xl)/(N-1);
x1 = x0:h:xl;
x2 = xl:h:xN;
x = [x1 x2];
t0 = 0;
T = 1;
%%%%%%%%%%%%%

%%%Constants%%%
b1 = 1;
b2 = 100;

%%% Initial condition %%%
x0 = -1/4;
r0 = 1/30;
u0 = @(x) exp(-(x0-x).^2/r0^2);
u0_t = 0;

%%% Exact solution %%%
%%%%%%%%%%%%%%%%%%%%%%

%%%% SBP4 %%%%%
[D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP4_D4(N, h);

%%%% Boundary conditions %%%%
%%%
L = [eN -e1; dN_1 -d1_1; b1*dN_2 -b2*d1_2; b1*dN_3 -b2*d1_3; e1 zeros(1,N); d1_2 zeros(1,N); zeros(1,N) eN; zeros(1,N) dN_2];


H = [H zeros(N); zeros(N) H];
HI = inv(H);
P = [eye(N) zeros(N); zeros(N) eye(N)] - HI*L'*inv(L*HI*L')*L;

A = (-P*[b1*D4 zeros(N); zeros(N) b2*D4]*P);

w0 = [u0(x1) u0(x2)];
w0_t = 0;

%%%% Solver %%%%
[w,k,t] = timestepper(t0,T,h,A,w0,w0_t);

%%%% Plot %%%%
figure(1);
for i = 1:20:length(t)
    plot(x1, w(1:N,i), '-b', x2, w(N+1:end,i), '-r');
%     hold on;
%     plot(x,u_exact(x,t(i)), 'r');
%     hold off;
    axis([-1 1 -1 1]);
    pause(0.00000001);
end

