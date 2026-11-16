function [Br, Bt, Bmag] = field_from_flux(fun, r, theta)
% Numerical B components from a flux function Xi(r,theta).
h = 1e-5;
dPsidr = (fun(r+h,theta) - fun(r-h,theta))/(2*h);
dPsidt = (fun(r,theta+h) - fun(r,theta-h))/(2*h);
Br = dPsidt/(r^2*sin(theta));
Bt = -dPsidr/(r*sin(theta));
Bmag = hypot(Br,Bt);
end
