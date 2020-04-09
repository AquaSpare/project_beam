%Wavefunction u_tt = c^2u_xx
pause on 
clear all
close all
%%%Domain%%%
N = [ 21 41 81 101 201 501 1001];
x0 = -1;
xl = 0;
xN = 1;
h = (xN-xl)./(N-1);
x1 = cell(length(N),1);
x2 = cell(length(N),1);
x = cell(length(N),1);
t0 = 0;
T = 0.4;

for i = 1:length(N)
    x1{i,1} = x0:h(i):xl;
    x2{i,1} = xl:h(i):xN;
    x{i,1} = [x1{i,1} x2{i,1}];
end
%%%%%%%%%%%%%

%%%Constant%%%

al = 1;
ar = 1;

bl = 1;
br = 8;

sl = sqrt(al*bl);
sr = sqrt(ar*br);

Trans = (2*sl)/(sl+sr);
R = (sl-sr)/(sl+sr);


cl = sqrt(bl/al);
cr = sqrt(br/ar);
% lambda = pi*c1/x2;
%%%%%%%%%%%

%%%Boundary conditions%%%
% u_x0 = 0; u_xN = 0; % Neumann

%%%Initial data%%%
% u0 = @(x) sin(lambda.*x/c) + 2*sin(2*lambda.*x/c) + sin(3*lambda.*x/c) + sin(4*lambda.*x/c) + sin(5*lambda.*x/c) + sin(6*lambda.*x/c);
xp = -1/4;
r0 = 1/30;
u0 = @(x) exp(-(xp-x).^2/r0^2);
u0_t = 0;

%%%Exact solution%%%
%u_exact = @(x,t) cos(lambda.*t)*sin(lambda.*x/c) + cos(2*lambda.*t)*2*sin(2*lambda.*x/c) + cos(3*lambda.*t)*sin(3*lambda.*x/c)+ cos(4*lambda.*t)*sin(4*lambda.*x/c) + cos(5*lambda.*t)*sin(5*lambda.*x/c) + cos(6*lambda.*t)*sin(6*lambda.*x/c);
D2 = cell(length(N),1);
H = cell(length(N),1);
HI = cell(length(N),1);
e1 = cell(length(N),1);
eN = cell(length(N),1);
d1 = cell(length(N),1);
dN = cell(length(N),1);

for i = 1:length(N)
    [D2{i,1}, H{i,1}, HI{i,1}, M{i,1}, e1{i,1}, eN{i,1}, d1{i,1}, dN{i,1}] = SBP4(N(i), h(i));
end

clt = cl;
crt = cr;

L = cell(length(N),1);
P = cell(length(N),1);
A = cell(length(N),1);

for i = 1:length(N)

    L{i,1} = [eN{i,1} -e1{i,1}; cl^2*dN{i,1} -crt^2*d1{i,1}; e1{i,1} zeros(1,N(i)); zeros(1,N(i)) eN{i,1}];
    H{i,1} = [H{i,1} zeros(N(i)); zeros(N(i)) H{i,1}];
    HI{i,1} = inv(H{i,1});

    P{i,1} = [eye(N(i)) zeros(N(i)); zeros(N(i)) eye(N(i))] - HI{i,1}*L{i,1}'*inv(L{i,1}*HI{i,1}*L{i,1}')*L{i,1};

    A{i,1} = (P{i,1}*[clt^2*D2{i,1} zeros(N(i)); zeros(N(i)) crt^2*D2{i,1}]*P{i,1});
end

w0 = cell(length(N),1);
w0_t = cell(length(N),1);
w = cell(length(N),1);
k = cell(length(N),1);
t = cell(length(N),1);
u_exact = cell(length(N),1);

for i = 1:length(N)
    w0{i,1} = [u0(x1{i,1}) u0(x2{i,1})];
    w0_t{i,1} = 0;
    [w{i,1},k{i,1},t{i,1}] = timestepper(t0,T,h(i),A{i,1},w0{i,1},w0_t{i,1});
end


for j = 1:length(N)
    u_exact{j,1} = zeros(length(w0{j,1}),length(t{j,1}));
    for i = 1:length(t{j,1})
        u_exact{j,1}(:,i) = exact_gausspulse(x{j,1},t{j,1}(i),xp,r0,cl,cr,R,Trans);
    end
end



error = zeros(length(N),1);

for i = 1:length(N)
    wtest = w{i,1}(:,end);
    utest = u_exact{i,1}(:,end);
    error(i) = (1./sqrt(N(i))).*norm(u_exact{i,1}(:,end)-w{i,1}(:,end));
end

loglog(h,error);
hold on 
xlabel('step size')
ylabel('error')
title('Convergence')