function w = homo_beam_ana_2(BC)
%HOMO_BEAM_ANA Summary of this function goes here
%   Detailed explanation goes here
if BC == 1
    %clamped
    Beta = 4.73004074486270402602404810;
    u = @(x) cosh(Beta.*x) - cos(Beta.*x) - (sinh(Beta)+sin(Beta))/(cosh(Beta)-cos(Beta))*(sinh(Beta.*x)-sin(Beta*x));
elseif BC == 2
    %free
     Beta = 4.73004074486270402602404810;
     u = @(x) cosh(Beta.*x) + cos(Beta.*x) - ((cos(Beta) -cosh(Beta))/(sin(Beta) - sinh(Beta)))*(sin(Beta.*x) +sinh(Beta.*x));
end

w = @(x,t) real(exp(-1i*Beta^2*t)*u(x));