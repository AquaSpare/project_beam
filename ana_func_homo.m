function y = ana_func_homo(x, BC)
b = x(1);
A1 = x(2);
A2 = x(3);
A3 = x(4);
A4 = x(5);

if(BC==1)%clamped
    y1 = A3 + A4;
    y2 = A1 + A2;
    y3 = A1*sinh(b) + A2*sin(b) + A3*cosh(b) + A4*cos(b);
    y4 = A1*cosh(b) + A2*cos(b) + A3*sinh(b) - A4*sin(b);
elseif(BC==2)%free
    y1 = A3 - A4;
    y2 = A1 - A2;
    y3 = A1*sinh(b) - A2*sin(b) + A3*cosh(b) - A4*cos(b);
    y4 = A1*cosh(b) - A2*cos(b) + A3*sinh(b) - A4*sin(b);
elseif(BC==3)%sliding
    y1 = A1 + A2;
    y2 = A1 - A2;
    y3 = A1*cosh(b) + A2*cos(b) + A3*sinh(b) - A4*sin(b);
    y4 = A1*cosh(b) - A2*cos(b) + A3*sinh(b) + A4*sin(b);
elseif(BC==4)%hinged
    y1 = A3 + A4;
    y2 = A3 - A4;
    y3 = A1*sinh(b) + A2*sin(b) + A3*cosh(b) + A4*cos(b);
    y4 = A1*sinh(b) - A2*sin(b) + A3*cosh(b) - A4*cos(b);
end
y = [y1;y2;y3;y4];