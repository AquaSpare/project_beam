%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2:a ordn.  SBP Finita differens         %%%
%%% operatorer framtagna av Ken Mattsson    %%%
%%% Har med 5te derivata                    %%%
%%% 2018 05 30                              %%% 
%%%                                         %%%
%%% 6 randpunkter, diagonal norm            %%%
%%%                                         %%%
%%% Datum: 2018-0-31             c          %%%
%%%                                         %%%
%%%                                         %%%
%%% H           (Normen)                    %%%
%%% D1          (approx f?rsta derivatan)   %%%
%%% D2          (approx andra derivatan)    %%%
%%% D3          (approx tredje derivatan)   %%%
%%% D4          (approx fj?rde derivatan)   %%%
%%% D5          (approx femte derivatan)    %%%
%%%                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% M?ste ange antal punkter (m) och stegl?ngd (h)
% Notera att dessa opetratorer ?r framtagna f?r att anv?ndas n?r
% vi har 3de och 4de derivator i v?r PDE
% I annat fall anv?nd de "traditionella" som har noggrannare
% randsplutningar f?r D1 och D2

% Vi b?rjar med normen. Notera att alla SBP operatorer delar samma norm,
% vilket ?r n?dv?ndigt f?r stabilitet 

H=diag(ones(m,1),0);H(1,1)=1/2;H(m,m)=1/2;


H=H*h;
HI=inv(H);


% First derivative SBP operator, 1st order accurate at first 6 boundary points

q1=1/2;
Q=q1*(diag(ones(m-1,1),1)-diag(ones(m-1,1),-1));

%Q=(-1/12*diag(ones(m-2,1),2)+8/12*diag(ones(m-1,1),1)-8/12*diag(ones(m-1,1),-1)+1/12*diag(ones(m-2,1),-2));


e_l=zeros(m,1);e_l(1)=1;
e_r=zeros(m,1);e_r(m)=1;


D1=HI*(Q-1/2*e_l*e_l'+1/2*e_r*e_r') ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Second derivative, 1st order accurate at first 6 boundary points
m1=-1;m0=2;
M=m1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+m0*diag(ones(m,1),0);M(1,1)=1;M(m,m)=1;
M=M/h;
    
d1_U=[-3/2 2 -1/2]/h;
d1_l=zeros(1,m);
d1_l(1:3)=d1_U;
d1_r=zeros(1,m);
d1_r(m-2:m)=fliplr(-d1_U);

D2=HI*(-M-e_l*d1_l+e_r*d1_r);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Third derivative, 1/h order accurate at first 6 boundary points

q2=1/2;q1=-1;
Q3=q2*(diag(ones(m-2,1),2)-diag(ones(m-2,1),-2))+q1*(diag(ones(m-1,1),1)-diag(ones(m-1,1),-1));   

%QQ3=(-1/8*diag(ones(m-3,1),3) + 1*diag(ones(m-2,1),2) - 13/8*diag(ones(m-1,1),1) +13/8*diag(ones(m-1,1),-1) -1*diag(ones(m-2,1),-2) + 1/8*diag(ones(m-3,1),-3));


Q3_U = [0 -0.13e2 / 0.16e2 0.7e1 / 0.8e1 -0.1e1 / 0.16e2; 0.13e2 / 0.16e2 0 -0.23e2 / 0.16e2 0.5e1 / 0.8e1; -0.7e1 / 0.8e1 0.23e2 / 0.16e2 0 -0.17e2 / 0.16e2; 0.1e1 / 0.16e2 -0.5e1 / 0.8e1 0.17e2 / 0.16e2 0;];
Q3(1:4,1:4)=Q3_U;
Q3(m-3:m,m-3:m)=flipud( fliplr( -Q3_U ) );
Q3=Q3/h^2;



d2_U=[1 -2 1;]/h^2;
d2_l=zeros(1,m);
d2_l(1:3)=d2_U;
d2_r=zeros(1,m);
d2_r(m-2:m)=fliplr(d2_U);

D3=HI*(Q3 - e_l*d2_l + e_r*d2_r +1/2*d1_l'*d1_l -1/2*d1_r'*d1_r ) ;

% Fourth derivative order?

m2=1;m1=-4;m0=6;
M4=m2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2))+m1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+m0*diag(ones(m,1),0);

%M4=(-1/6*(diag(ones(m-3,1),3)+diag(ones(m-3,1),-3) ) + 2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2)) -13/2*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1)) + 28/3*diag(ones(m,1),0));

M4_U=[0.13e2 / 0.10e2 -0.12e2 / 0.5e1 0.9e1 / 0.10e2 0.1e1 / 0.5e1; -0.12e2 / 0.5e1 0.26e2 / 0.5e1 -0.16e2 / 0.5e1 0.2e1 / 0.5e1; 0.9e1 / 0.10e2 -0.16e2 / 0.5e1 0.47e2 / 0.10e2 -0.17e2 / 0.5e1; 0.1e1 / 0.5e1 0.2e1 / 0.5e1 -0.17e2 / 0.5e1 0.29e2 / 0.5e1;];


M4(1:4,1:4)=M4_U;

M4(m-3:m,m-3:m)=flipud( fliplr( M4_U ) );
M4=M4/h^3;
    
d3_U=[-1 3 -3 1;]/h^3;
d3_l=zeros(1,m);
d3_l(1:4)=d3_U;
d3_r=zeros(1,m);
d3_r(m-3:m)=fliplr(-d3_U);

D4=HI*(M4-e_l*d3_l+e_r*d3_r  + d1_l'*d2_l-d1_r'*d2_r);


% Fifth derivative, 1/h^2 order accurate at first 3 boundary points

q3=1/2;q2=-2;q1=5/2;
Q5=q3*(diag(ones(m-3,1),3)-diag(ones(m-3,1),-3))+q2*(diag(ones(m-2,1),2)-diag(ones(m-2,1),-2))+q1*(diag(ones(m-1,1),1)-diag(ones(m-1,1),-1));   

%QQ3=(-1/8*diag(ones(m-3,1),3) + 1*diag(ones(m-2,1),2) - 13/8*diag(ones(m-1,1),1) +13/8*diag(ones(m-1,1),-1) -1*diag(ones(m-2,1),-2) + 1/8*diag(ones(m-3,1),-3));


Q5_U = [0 0.1e1 / 0.2e1 -1 0.1e1 / 0.2e1; -0.1e1 / 0.2e1 0 2 -2; 1 -2 0 0.5e1 / 0.2e1; -0.1e1 / 0.2e1 2 -0.5e1 / 0.2e1 0;];

Q5(1:4,1:4)=Q5_U;
Q5(m-3:m,m-3:m)=flipud( fliplr( -Q5_U ) );
Q5=Q5/h^4;

d4_U=[1 -4 6 -4 1;]/h^4;
d4_l=zeros(1,m);
d4_l(1:5)=d4_U;
d4_r=zeros(1,m);
d4_r(m-4:m)=fliplr(d4_U);

D5=HI*(Q5 - e_l*d4_l + e_r*d4_r + d1_l'*d3_l -d1_r'*d3_r -1/2*d2_l'*d2_l + 1/2*d2_r'*d2_r) ;


%L=h*(m-1);

% x1=linspace(0,L,m)';
% x2=x1.^2/factorial(2);
% x3=x1.^3/factorial(3);
% x4=x1.^4/factorial(4);
% x5=x1.^5/factorial(5);
% x6=x1.^6/factorial(6);
% x7=x1.^7/factorial(7);
% x8=x1.^8/factorial(8);
% x9=x1.^9/factorial(9);
% 
% 
% x0=x1.^0;

