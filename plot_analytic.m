clear;
close all;
BC = 1;
a1 = 1;
a2 = 100;

timestep = 0.0005;
T = 2;
t = 0:timestep:T;

x0 = -1; xl = 0; xN = 1;
N = 81;
x1 = x0:(xl-x0)/(N-1):xl;
x2 = xl:(xN-xl)/(N-1):xN;

[u,v] = beam_ana(BC, a1, a2);
[u2,v2] = beam_ana_2(BC, a1, a2);
ymax = max(abs(u2(x1,0)));

for i = 1:length(t)
    plot(x1,u(x1,t(i)), 'b', x2, v(x2,t(i)), 'r');
    hold on;
    plot(x1,u2(x1,t(i)), 'g', x2, v2(x2,t(i)), 'y');
    hold off;
    axis([x0 xN -ymax ymax]);
    pause(0.00000001);
end