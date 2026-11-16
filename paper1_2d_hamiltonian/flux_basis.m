function [Phi, dPhidr, dPhidt] = flux_basis(r, theta, n)
% Axisymmetric internal flux basis (zonal harmonic):
%   Phi_n(r,theta) = - sin(theta) P_n^1(cos(theta)) / (n r^n)
% For n=1 this is Phi_1 = sin^2(theta)/r, i.e. the dipole flux function
% divided by g10.
%
% Units: r in R_E, Phi in R_E (with R_E=1 normalization), angles in rad.
x = cos(theta);
P = legendre(n, x);
Pn  = P(1);
Pn1 = P(2);
s = sin(theta);
Phi = -(s .* Pn1) ./ (n .* r.^n);
dPhidr  = -n .* Phi ./ r;
dPhidt  = (n + 1) .* s .* Pn ./ r.^n;
end
