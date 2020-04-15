%exempel
%beam_homogenous(81,0,1,0.18,1,4,4,500,0.0001,1)

function [w,t,x,h,k] = beam_homogenous(N,x0,xl,T,b,order,BC,plotspeed,k_ratio,plotornot, timestepperversion)
close all
pause on

%mesh and domain
h = (xl-x0)/(N-1);
x = x0:h:xl;

%xp = 0.5;
%r0 = 0.05;


%initial condition
if b == 1 && x0 == 0 && xl == 1
    Beta = 4.73004074486270402602404810;
    u0 = @(x) cosh(Beta.*x) + cos(Beta.*x) - ((cos(Beta) -cosh(Beta))/(sin(Beta) - sinh(Beta)))*(sin(Beta.*x) +sinh(Beta.*x));
    w0_t = 0;
    w0 = u0(x);
    special_case = 1;
else
    xp = 0.5;
    r0 = 0.05;
    u0 = @(x) exp(-(xp-x).^2/r0^2);
    w0_t = 0;
    w0 = u0(x);
    special_case = 0;
end

if order == 2
    [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP2_D4(N, h);
elseif order == 4
    [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP4_D4(N, h);
elseif order == 6
    [D4, H, HI, M, e1, eN, d1_1, dN_1, d1_2, dN_2, d1_3, dN_3] = SBP6_D4(N, h);
end

%HI = inv(H);

%------------Boundary conditions--------------
if BC == 1
    %------Clamped-------
    L = [e1;eN;d1_1;dN_1];
elseif BC == 2
    %------Hinged--------
    L = [e1;eN;d1_2;dN_2];
elseif BC == 3
    %------Sliding-------
    L = [d1_1;dN_1;d1_3;dN_3];
elseif BC == 4
    %-------Free---------
    L = [d1_2;dN_2;d1_3;dN_3];
end
    
    
P = eye(N)-HI*L'*inv(L*HI*L')*L;
A = -b.*P*D4*P;

if timestepperversion == 1
    [w,k,t] = timestepper(T,h,A,w0,w0_t,k_ratio);
elseif timestepperversion == 2
    [w,k,t] = timestepperv2(T,h,A,w0,w0_t,k_ratio);
end
%exact soliution for special case
u_exact = @(x,t) real(exp(-1i*(Beta^2)*t)*u0(x));

if plotornot == 1 && special_case == 0
    for i = 1:plotspeed:length(t)
        plot(x,w(:,i))
        axis([0 1 -5 5])
        pause(0.00000001)
    end
elseif plotornot == 1 && special_case == 1
    for i = 1:plotspeed:length(t)
        plot(x,w(:,i),'*')
        hold on
        plot(x,u_exact(x,t(i)))
        axis([0 1 -5 5])
        pause(0.00000001)   
        hold off
    end
end
end

