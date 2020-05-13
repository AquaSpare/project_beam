clear;
close all;
BC = 1;
a1 = 1;
a2 = 100;

timestep = 0.0005;
T = 0.09;
t = 0:timestep:T;

x0 = -1; xl = 0; xN = 1;
N = 81;
x1 = x0:(xl-x0)/(N-1):xl;
x2 = xl:(xN-xl)/(N-1):xN;

[u,v] = beam_ana(BC, a1, a2);
ymax = max(abs([u(x1,0) v(x2,0)]));

for i = 1:length(t)
    plot(x1,u(x1,t(i)), 'b', x2, v(x2,t(i)), 'r');
    axis([x0 xN -ymax ymax]);
    xlabel('x','FontSize',16);
    title('Analytic solution, clamped','FontSize',18);
    pause(0.00000001);
end
