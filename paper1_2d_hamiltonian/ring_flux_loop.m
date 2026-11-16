function [Psi, Aphi] = ring_flux_loop(R, Z, a, Bc)
% Flux function Psi = R*Aphi for a circular current loop of radius a (R_E)
% in the equatorial plane.  Bc is the field at the loop centre (nT).
% The field of the loop is B_R = -(1/R) dPsi/dZ, B_Z = (1/R) dPsi/dR.
if R < 1e-10
    Psi = 0; Aphi = 0; return;
end
beta2 = (a + R)^2 + Z^2;
k2 = 4*a*R / beta2;
k2 = min(max(k2, 0), 1 - 1e-12);
[K,E] = ellipke(k2);
f = (4*a/sqrt(beta2)) * (((2 - k2)*K - 2*E)/k2);
C = Bc * a / (2*pi);
Aphi = C * f;
Psi = R * Aphi;
end
