%exempel
%beam_homogenous(81,0,1,0.18,1,4,4,500,0.0001,1,1,1);

function [w,t,x,h,k,error] = beam_homogenous(N,x0,xN,T,a,order,BC,plotspeed,k_ratio,IC,plotornot, timestepperversion, timestepperOrder)
close all
pause on

%mesh and domain
h = (xN-x0)/(N-1);
x = x0:h:xN;


%%% Initial condition %%%
if(IC == 0)
    xr = -1/4;
    r0 = 1/30;
    u = @(x) exp(-(xr-x).^2/r0^2);
    u0 = u(x);
else 
    u = homo_beam_ana_2(a, BC);
    u0 = u(x,0);
    u0_t = zeros(1,length(u0));
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
    %-------Free---------
    L = [d1_2;dN_2;d1_3;dN_3];
elseif BC == 3
    %------Sliding-------
    L = [d1_1;dN_1;d1_3;dN_3];
elseif BC == 4
    %------Hinged--------
    L = [e1;eN;d1_2;dN_2];
end
    
    
P = eye(N)-HI*L'*inv(L*HI*L')*L;
A = -a.*P*D4*P;

if timestepperversion == 1
    [w,k,t] = timestepper(T,h,A,u0,u0_t,k_ratio,timestepperOrder);
    error = 0;
elseif timestepperversion == 2
    [w,k,t] = timestepperv2(T,h,A,u0,u0_t,k_ratio,timestepperOrder);
    error = 0;
elseif timestepperversion == 3
    [w,k,t,error] = timestepperv3_enBalk(T,h,A,u0,u0_t,k_ratio,BC,a,x0,xN);
end

%%% Plot
if plotornot == 1
    figure(1);
    if(IC ~= 0)
        ymax = max(abs(u(x,0)));
    end
    for i = 1:plotspeed:length(t)
        if(IC ~= 0)
            plot(x, w(1:N,i), '*b');
            hold on;
            plot(x,u(x,t(i)), 'r');
            axis([x0 xN -ymax ymax]);
            hold off;
        else
            plot(x1, w(1:N,i), 'b');
            axis([x0 x1 -1 1]);
        end
        pause(0.00000001);
    end
end
end

