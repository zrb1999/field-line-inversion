function [Phi, dPhidr, dPhidt] = flux_basis_ext(r, theta, n)
% Axisymmetric EXTERNAL flux basis (positive power of r):
%   Phi_n^ext(r,theta) = - sin(theta) P_n^1(cos(theta)) r^(n+1) / (n+1)
% For n=1 this is Phi_1^ext = r^2 sin^2(theta)/2, i.e. a uniform field
% along the symmetry axis (B_z = G_1).
%
% Units: r in R_E, Phi in R_E^2 (with R_E=1 normalization), angles in rad.
x = cos(theta);
P = legendre(n, x);
Pn  = P(1);
Pn1 = P(2);
s = sin(theta);
Phi = -(s .* Pn1) .* r.^(n+1) ./ (n+1);
dPhidr = (n+1) .* Phi ./ r;
dPhidt = n .* r.^(n+1) .* s .* Pn;
end
