function y = func(x,a1,a2,BC)
b1 = x(1);
b2 = x(2);
A1 = x(3);
A2 = x(4);
A3 = x(5);
A4 = x(6);
A5 = x(7);
A6 = x(8);
A7 = x(9);
A8 = x(10);

if(BC == 1) % Clamped
    y1 = A1*sinh(-b1) + A2*sin(-b1) + A3*cosh(-b1) + A4*cos(-b1);
    y2 = A1*cosh(-b1) + A2*cos(-b1) + A3*sinh(-b1) - A4*sin(-b1);
    y3 = A5*sinh(b2)  + A6*sin(b2)  + A7*cosh(b2)  + A8*cos(b2);
    y4 = A5*cosh(b2)  + A6*cos(b2)  + A7*sinh(b2)  - A8*sin(b2);
elseif(BC == 2) % Free
    y1 = A1*sinh(-b1) - A2*sin(-b1) + A3*cosh(-b1) - A4*cos(-b1);
    y2 = A1*cosh(-b1) - A2*cos(-b1) + A3*sinh(-b1) + A4*sin(-b1);
    y3 = A5*sinh(b2)  - A6*sin(b2)  + A7*cosh(b2)  - A8*cos(b2);
    y4 = A5*cosh(b2)  - A6*cos(b2)  + A7*sinh(b2)  + A8*sin(b2);
elseif(BC == 3) % Sliding
    y1 = A1*cosh(-b1) + A2*cos(-b1) + A3*sinh(-b1) - A4*sin(-b1);
    y2 = A1*cosh(-b1) - A2*cos(-b1) + A3*sinh(-b1) + A4*sin(-b1);
    y3 = A5*cosh(b2)  + A6*cos(b2)  + A7*sinh(b2)  - A8*sin(b2);
    y4 = A5*cosh(b2)  - A6*cos(b2)  + A7*sinh(b2)  + A8*sin(b2);
elseif(BC == 4)% Hinged
    y1 = A1*sinh(-b1) + A2*sin(-b1) + A3*cosh(-b1) + A4*cos(-b1);
    y2 = A1*sinh(-b1) - A2*sin(-b1) + A3*cosh(-b1) - A4*cos(-b1);
    y3 = A5*sinh(b2)  + A6*sin(b2)  + A7*cosh(b2)  + A8*cos(b2);
    y4 = A5*sinh(b2)  - A6*sin(b2)  + A7*cosh(b2)  - A8*cos(b2);
end

y5 = sqrt(a1)*b1^2 - sqrt(a2)*b2^2;
y6 = A3 + A4 - (A7 + A8);
y7 = b1*(A1 + A2) - b2*(A5 + A6);
y8 = a1*b1^2*(A3 - A4) - a2*b2^2*(A7 - A8);
y9 = a1*b1^3*(A1 - A2) - a2*b2^3*(A5 - A6);

y = [y1;y2;y3;y4;y5;y6;y7;y8;y9];
