close all;
N = [21 41 81 161 321];
x0 = 0;
xN = 1;
h = (xN-x0)./(N-1);
k = zeros(1,length(N));
k2 = zeros(1,length(N));
t0 = 0;
T = 0.7;
c = 1;
lambda = pi*c/xN;

%%%%%% Initial data %%%%%%
u_x0 = 0; u_xN = 0; % Neumann
u0 = @(x) cos(lambda.*x/c) + 2*cos(2*lambda.*x/c) + cos(3*lambda.*x/c) + cos(4*lambda.*x/c) + cos(5*lambda.*x/c) + cos(6*lambda.*x/c);
u0_t = 0;
u_exact = @(x,t) cos(lambda.*t)*cos(lambda.*x/c) + cos(2*lambda.*t)*2*cos(2*lambda.*x/c) + cos(3*lambda.*t)*cos(3*lambda.*x/c)+ cos(4*lambda.*t)*cos(4*lambda.*x/c) + cos(5*lambda.*t)*cos(5*lambda.*x/c) + cos(6*lambda.*t)*cos(6*lambda.*x/c);
%%%%%%%%%%%%%%%%%%%%%%%%%%

e1 = cell(1,length(N));
eN = cell(1,length(N));
d1 = cell(1,length(N));
dN = cell(1,length(N));
H = cell(1,length(N));
M = cell(1,length(N));
D2 = cell(1,length(N));
u = cell(1, length(N));
error = cell(1, length(N));
x = cell(1, length(N));
t = cell(1, length(N));
for j = 1:length(N)
   e1{1,j} = zeros(1,N(j));
   eN{1,j} = zeros(1,N(j));
   d1{1,j} = zeros(1,N(j));
   dN{1,j} = zeros(1,N(j));
   e1{1,j}(1,1) = 1;   
   eN{1,j}(1,end) = 1;
   %%%%%% second-order D2%%%%%% 
%    d1{1,j}(1,1:3) = 1/h(j)*[-3/2 2 -1/2];
%    dN{1,j}(1,(end-2):end) = 1/h(j)*[1/2 -2 3/2];
%    H{1,j} = h(j)*diag(ones(1,N(j)));
%    H{1,j}(1,1) = H{1,j}(1,1)/2; 
%    H{1,j}(end,end) = H{1,j}(end,end)/2;
%    M{1,j} = (1/h(j))*(zeros(N(j))  + diag(2*ones(N(j) ,1)) + diag(-ones(N(j)-1,1), 1) + diag(-ones(N(j)-1,1), -1));  
%    M{1,j}(1,1) = M{1,j}(1,1)/2; 
%    M{1,j}(end,end) = M{1,j}(end,end)/2;

   %%%%%% fourth-order D2%%%%%%
   d1{1,j}(1,1:4) = 1/h(j)*[11/6 -3 3/2 -1/3];
   dN{1,j}(1,(end-3):end) = 1/h(j)*[1/3 -3/2 3 -11/6];
   H{1,j} = diag(ones(1,N(j)));
   H{1,j}(1,1) = 17/48; H{1,j}(end,end) = 17/48;
   H{1,j}(2,2) = 59/48; H{1,j}(end-1,end-1) = 59/48;
   H{1,j}(3,3) = 43/48; H{1,j}(end-2,end-2) = 43/48;
   H{1,j}(4,4) = 49/48; H{1,j}(end-3,end-3) = 49/48;
   H{1,j} = h(j) * H{1,j};
   M{1,j} = (zeros(N(j))  + diag(5/2*ones(N(j) ,1)) + diag(-4/3*ones(N(j)-1,1), 1) + diag(-4/3*ones(N(j)-1,1), -1) + diag(1/12*ones(N(j)-2,1), 2) + diag(1/12*ones(N(j)-2,1), -2));  
   M_U1 = [9/8 -59/48 1/12 1/48]; M_U2 = [-59/48 59/24 -59/48 0]; M_U3 = [1/12 -59/48 55/24 -59/48]; M_U4 = [1/48 0 -59/48 59/24];
   M{1,j}(1,1:4) = M_U1; M{1,j}(2,1:4) = M_U2; M{1,j}(3,1:4) = M_U3; M{1,j}(4,1:4) = M_U4;
   M{1,j}(end,(end-3):end) = fliplr(M_U1); M{1,j}(end-1,(end-3):end) = fliplr(M_U2); M{1,j}(end-2,(end-3):end) = fliplr(M_U3); M{1,j}(end-3,(end-3):end) = fliplr(M_U4);
   M{1,j} = 1/h(j) * M{1,j};
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   D2{1,j} = H{1,j}\(-M{1,j}-e1{1,j}'*d1{1,j}+eN{1,j}'*dN{1,j});
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%%%%%% timestep %%%%%%%
%    k(j) = 2*h(j)/sqrt(max(abs(eig(c^2.*(D2{1,j} + H{1,j}\(e1{1,j}'*d1{1,j}) - H{1,j}\(eN{1,j}'*dN{1,j})))))); 
   k(j) = 0.01*h(j)/c; 
   %%%%%%%%%%%%%%%%%%%%%%%%
   
   x{1,j} = linspace(x0, xN, N(j));
   t{1,j} = linspace(t0, T, (T-t0)/k(j));
   u{1,j} = zeros(length(x{1,j}), length(t{1,j}));
   error{1,j} = zeros(1, length(t{1,j}));
   %%%%% Initial data %%%%%%%%%
   u{1,j}(:,1) = u0(x{1,j});
   u{1,j}(:,2) = k(j)*u0_t + k(j)^2/2*c^2*(D2{1,j}*u{1,j}(:,1)+ H{1,j}\(e1{1,j}'*(d1{1,j}*u{1,j}(:,1)-u_x0)) - H{1,j}\(eN{1,j}'*(dN{1,j}*u{1,j}(:,1)-u_xN))) + u{1,j}(:,1);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pause on;
   for n = 3:length(t{1,j})
       u{1,j}(:,n) = k(j)^2*c^2.*(D2{1,j}*u{1,j}(:,n-1) + H{1,j}\(e1{1,j}'*(d1{1,j}*u{1,j}(:,n-1)-u_x0)) - H{1,j}\(eN{1,j}'*(dN{1,j}*u{1,j}(:,n-1)-u_xN))) + 2*u{1,j}(:,n-1) - u{1,j}(:,n-2);
       error{1,j}(1,n) = norm(-u{1,j}(:,n)'+ u_exact(x{1,j},t{1,j}(n)),2)/sqrt(N(j));
%        figure(j);
%        plot(x{1,j},u{1,j}(:,n),'*b');
%        hold on;
%        plot(x{1,j},u_exact(x{1,j},t{1,j}(n)),'r');
%        axis([x0 xN -4 7]);
%        pause(0.0001);
%        hold off;
   end   
   figure(j);
   plot(x{1,j},u{1,j}(:,end),'*b');
   hold on;
   plot(x{1,j},u_exact(x{1,j},t{1,j}(end)),'r');
   hold off;
end
figure(j+1);
err = zeros(1, length(N));
for j=1:length(err)
   err(j) = norm(-u{1,j}(:,end)'+ u_exact(x{1,j},t{1,j}(end)),2)/sqrt(N(j));
end
p = polyfit(log(h),log(err),1);
y = polyval(p,log(h));
loglog(h, exp(y),'r',h,err, 'b*');
legend([num2str(p(1)),'log(h) + ', num2str(p(2))],'Error');
xlabel('log(h)');
ylabel('log(error)');
title('Error estimation');
