clear all
close all
pause on

T = 50;
a = 1;

N = 81;
order = 4;
kratio = 0.1;
timestepperOrder = 4;

BC = 3;
x0 = 0;
xN = 1;

u_exact = homo_beam_ana_2(a, BC);
error = cell(2,1);
u = cell(1,2);

[u{1,1}, t, x, h, k, error{1,1}] = SAT_enBalk(N,x0,xN,T,a,order, BC, 500, kratio, 1, 0, 3,timestepperOrder);
[u{1,2},t,x,h,k, error{2,1}] = beam_homogenous(N,x0,xN,T,a,order,BC,500,kratio,1, 0, 3,timestepperOrder);

figure(1);
plot(t, error{1,1});
hold on
plot(t, error{2,1});
legend('SAT', 'Proj');
hold off

figure(2);
plot(x, u{1,1}(:,3), 'b*');
hold on;
plot(x0:h:xN, u_exact(x0:h:xN, T), 'r');