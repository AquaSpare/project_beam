clear all
close all

c = 1;

xr = 1;
xl = 0;
L = xr-xl;

N = [21 41 81 161 321];
h = (xr-xl)./(N-1);
%tummregel f�r k
k = h./(2*c);
J = length(N);

H = cell(J,1);
M = cell(J,1);
d1 = cell(J,1);
dm = cell(J,1);
e1 = cell(J,1);
em = cell(J,1);

for i = 1:J
%---------------------- SECOND ORDER -----------------------%

H{i,1} = diag(ones(N(i),1));
H{i,1}(1,1) = 1/2;
H{i,1}(end,end) = 1/2;
H{i,1} = h(i).*H{i,1};

M{i,1} = diag(2.*ones(N(i),1)) + diag(-1.*ones(N(i)-1,1),1) + diag(-1.*ones(N(i)-1,1),-1);
M{i,1}(1,1) = 1;
M{i,1}(end,end) = 1;
M{i,1} = h(i).*M{i,1};

e1{i,1} = zeros(1,J(i));
e1{i,1}(1,1) = 1;
em{i,1} = zeros(1,J(i));
em{i,1}(1,end) = 1;

d1{i,1} = zeros(1,J(i));
d1{i,1}(1,1:3) = [-3/2 2 1/2];
d1{i,1} = d1{i,1}.*(1/h(i));
dm{i,1} = zeros(1,J(i));
dm{i,1}(1,(end-2):end) = [1/2 -2 3/2];
dm{i,1} = dm{i,1}.*(1/h(i));
end