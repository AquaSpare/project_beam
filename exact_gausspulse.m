
function u = exact_gausspulse(x,t,xp,r0,cl,cr,R,T)

t_star = -xp/cl;

if t<t_star
    u = 0.5*exp(-(xp-(x+cl*t))^2/(r0^2)) + 0.5*exp(-(xp-(x-cl*t))^2/(r0^2)); 
else
    u = 
end

        