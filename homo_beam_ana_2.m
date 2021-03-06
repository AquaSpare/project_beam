function w = homo_beam_ana_2(a, BC)
%HOMO_BEAM_ANA Summary of this function goes here
%   Detailed explanation goes here

if BC == 1%clamped
    % a = 1
    Beta = 4.73004074486270402602404810;
    u = @(x) cosh(Beta.*x) - cos(Beta.*x) - (sinh(Beta)+sin(Beta))/(cosh(Beta)-cos(Beta))*(sinh(Beta.*x)-sin(Beta*x));
    w = @(x,t) real(exp(-1i*Beta^2*t)*u(x));
elseif BC == 2%free
    % a = 1
     Beta = 4.73004074486270402602404810;
     u = @(x) cosh(Beta.*x) + cos(Beta.*x) - ((cos(Beta) -cosh(Beta))/(sin(Beta) - sinh(Beta)))*(sin(Beta.*x) +sinh(Beta.*x));
     w = @(x,t) real(exp(-1i*Beta^2*t)*u(x));
elseif(BC == 3)%sliding
    n = 2;
    b = pi*n;
    w = @(x,t) cos(sqrt(a)*b^2*t)*cos(b*x);
elseif(BC == 4)%hinged
    n = 2;
    b = pi*n;
    w = @(x,t) cos(sqrt(a)*b^2*t)*sin(b*x);
end

