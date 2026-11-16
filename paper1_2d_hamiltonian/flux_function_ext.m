function [Xi, Br, Bt, Bmag] = flux_function_ext(r, theta, gint, Gext)
% Flux function and field for internal + external zonal harmonics:
%   Xi = sum_n gint(n) Phi_n^int + sum_n Gext(n) Phi_n^ext
%   B_r = (1/(r^2 sin theta)) dXi/dtheta
%   B_theta = -(1/(r sin theta)) dXi/dr
Ni = numel(gint);
Ne = numel(Gext);
Xi = 0; dXi_dr = 0; dXi_dt = 0;

for n = 1:Ni
    [Phi, dPr, dPt] = flux_basis(r, theta, n);
    Xi     = Xi     + gint(n) * Phi;
    dXi_dr = dXi_dr + gint(n) * dPr;
    dXi_dt = dXi_dt + gint(n) * dPt;
end
for n = 1:Ne
    [Phi, dPr, dPt] = flux_basis_ext(r, theta, n);
    Xi     = Xi     + Gext(n) * Phi;
    dXi_dr = dXi_dr + Gext(n) * dPr;
    dXi_dt = dXi_dt + Gext(n) * dPt;
end

Br = dXi_dt ./ (r.^2 .* sin(theta));
Bt = -dXi_dr ./ (r .* sin(theta));
Bmag = hypot(Br, Bt);
end
