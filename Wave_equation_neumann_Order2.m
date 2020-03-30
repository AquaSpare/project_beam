clear all
close all
pause on
N = [321 641 1281 2561];
 
%u_tt = c^2u_xx
x0 = 0;
xN = 1;
h = (xN-x0)./(N-1);
k = zeros(1,length(N));
k2 = zeros(1,length(N));
k3 = zeros(1,length(N));
t0 = 0;
T = 2;
c = 1*[1 1 1 1 1 1 1];
% c = h./4
lambda = pi.*c./xN;
%%%%%% Initial data %%%%%%
% Homogen Neumann
u_x0 = 1;
u_xN= 0;
%different solutions
% u0 = @(x,lambda,c) cos(lambda*x/c);
% u0_t = 0;
% u_exact = @(x,t,lambda,c) cos(lambda.*t).*cos(lambda*x/c);
% %%%%%%
u0 = @(x,lambda,c) cos(lambda.*x/c) + 2*cos(2*lambda.*x/c) + cos(3*lambda.*x/c) + cos(4*lambda.*x/c) + cos(5*lambda.*x/c) + cos(6*lambda.*x/c);
u0_t = 0;
u_exact = @(x,t,lambda,c) cos(lambda.*t)*cos(lambda.*x/c) + cos(2*lambda.*t)*2*cos(2*lambda.*x/c) + cos(3*lambda.*t)*cos(3*lambda.*x/c)+ cos(4*lambda.*t)*cos(4*lambda.*x/c) + cos(5*lambda.*t)*cos(5*lambda.*x/c) + cos(6*lambda.*t)*cos(6*lambda.*x/c);
%%%%%%%%%%%%%%%%%%%%%%%%%%

e1 = cell(1,length(N));
eN = cell(1,length(N));
d1 = cell(1,length(N));
dN = cell(1,length(N));
H = cell(1,length(N));
M = cell(1,length(N));
D2 = cell(1,length(N));
u = cell(1, length(N));
x = cell(1, length(N));
t = cell(1, length(N));
for j = 1:length(N)
   %%%%%% second-order D2%%%%%%
   e1{1,j} = zeros(1,N(j));
   eN{1,j} = zeros(1,N(j));
   d1{1,j} = zeros(1,N(j));
   dN{1,j} = zeros(1,N(j));
   e1{1,j}(1,1) = 1;   
   eN{1,j}(1,end) = 1;
   d1{1,j}(1,1:4) = 1/h(j)*[11/6 -3 3/2 -1/3];
   dN{1,j}(1,(end-3):end) = 1/h(j)*-[-1/3 3/2 -3 11/6];
   %%%H-Matris%%%
   H{1,j} = h(j).*diag(ones(1,N(j)));
   H{1,j}(1,1) = H{1,j}(1,1)*(17/48); H{1,j}(end,end) = H{1,j}(1,1);
   H{1,j}(2,2) = H{1,j}(2,2)*(59/48); H{1,j}(end-1,end-1) = H{1,j}(2,2);
   H{1,j}(3,3) = H{1,j}(3,3)*(43/48); H{1,j}(end-2,end-2) = H{1,j}(3,3);
   H{1,j}(4,4) = H{1,j}(4,4)*(49/48); H{1,j}(end-3,end-3) = H{1,j}(4,4);
   
   %%%M-Matris%%%
   M{1,j} = (zeros(N(j))  + diag((5/2)*ones(N(j) ,1)) + diag(-(4/3)*ones(N(j)-1,1), 1) + diag(-(4/3)*ones(N(j)-1,1), -1)+ diag((1/12)*ones(N(j)-2,1), 2) + diag((1/12)*ones(N(j)-2,1), -2));  
   M_U = [9/8 -59/48 1/12 1/48];
   M{1,j}(1,1:4) = M_U; M{1,j}(end,end-3:end)=rot90(M_U, 2);
   M{1,j}(1:4,1) = rot90(M_U,3); M{1,j}(end-3:end,end) = rot90(M_U,1);
   M{1,j}(2,2) = 59/24; M{1,j}(3,3)= 55/24; M{1,j}(4,4) =59/24;
   M{1,j}(2,3) = -59/48; M{1,j}(3,2) = -59/48; M{1,j}(4,3)=-59/48; M{1,j}(3,4) = -59/48;
   M{1,j}(end-1,end-2) = -59/48; M{1,j}(end-2,end-1)=-59/48; M{1,j}(end-3,end-2)=-59/48; M{1,j}(end-2,end-3)=-59/48;
   for i =1:3
       M{1,j}(end-i,end-i)=M{1,j}(i+1,i+1);
   end
   M{1,j}(4,2) = 0; M{1,j}(2,4) = 0; M{1,j}(end-3,end-1)=0; M{1,j}(end-1,end-3) = 0;
   M{1,j} = (1/h(j)).*M{1,j};
   
   %%%second derivative%%%
   D2{1,j} = H{1,j}\(-M{1,j}-e1{1,j}'*d1{1,j}+eN{1,j}'*dN{1,j});
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%%%%%% timestep %%%%%%%
   k(j) = (2*h(j)/(sqrt(max(abs(eig(c(j)^2*D2{1,j}*h(j)^2))))));
   %k2(j) = h(j)/c;
%    k(j) = 0.00001
   %%%%%%%%%%%%%%%%%%%%%%%%
   
   x{1,j} = x0:h(j):xN;
   t{1,j} = t0:k(j):T;
   u{1,j} = zeros(length(x{1,j}), length(t{1,j}));
   %%%%% Initial data %%%%%%%%%
   u{1,j}(:,1) = u0(x{1,j},lambda(j),c(j));
   u{1,j}(:,2) = k(j)*u0_t + ((k(j)^2)/2)*c(j)^2*D2{1,j}*u{1,j}(:,1) +u{1,j}(:,1) ;
 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   shg
   figure(10)
   for n = 3:length(t{1,j})
%        u{1,j}(:,n) = k(j)^2*c(j)^2.*D2{1,j}*u{1,j}(:,n-1) +2*u{1,j}(:,n-1) -u{1,j}(:,n-2) ;
       u{1,j}(:,n) = k(j)^2.*c(j)^2.*(D2{1,j}*u{1,j}(:,n-1) + H{1,j}*(e1{1,j}'*(d1{1,j}*u{1,j}(:,n-1)-u_x0)) - H{1,j}*(eN{1,j}'*(dN{1,j}*u{1,j}(:,n-1)-u_xN))) + 2*u{1,j}(:,n-1) - u{1,j}(:,n-2);
        %%%plottning%%%
       plot(x{1,j},u{1,j}(:,n),'*b');
       axis([x0 xN -1 1]);
       hold on
       plot(x{1,j},u_exact(x{1,j},t{1,j}(n),lambda(j),c(j)),'r');
       axis([x0 xN -7 7]);
       hold off
       pause(0.001)

       %u(n+1) = k^2*c^2*D2*u(n)+ 2*u(n) - u(n-1);
       %u_tt = 1/k^2*(u(n+1) - 2*u(n) + u(n-1)) = c^2*D2*u(n) = c^2*u_xx;
   end
   figure(j)
   hold on
   plot(x{1,j},u{1,j}(:,end),'b')
   plot(x{1,j},u_exact(x{1,j},t{1,j}(end),lambda(j),c(j)),'r')
   hold off

   error(j) = (1/sqrt(N(j))).*norm(u_exact(x{1,j},t{1,j}(end),lambda(j),c(j))'-u{1,j}(:,end),2);
end
p = polyfit(log(h),log(error),1);
y = polyval(p,log(h));
loglog(h, exp(y),'r',h,error, 'b*');
legend([num2str(p(1)),'log(h) + ', num2str(p(2))],'Error');
xlabel('log(h)');
ylabel('log(error)');
title('Error estimation');