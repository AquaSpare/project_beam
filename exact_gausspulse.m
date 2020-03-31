
function u = exact_gausspulse(x,t,xp,r0,cl,cr,R,T)

t_star = -xp/cl;

if t<t_star
    u = 0.5*exp(-((xp-(x+cl*t)).^2)./(r0^2)) + 0.5*exp(-((xp-(x-cl*t)).^2)./(r0^2)); 
else
    xl = x(x<0);
    xr = x(x>= 0);
    
    ul = 0.5*exp(-((xp-(xl+cl*t)).^2)./(r0^2)) + (R/2)*exp(-((xl+cl*(t-t_star)).^2)./(r0^2));
    ur = (T/2)*exp(-((cl*(xr-(cr*(t-t_star)))).^2)./(cr*r0)^2);
    
    u = [ul ur];
end

        