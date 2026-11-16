function Xi = flux_mixed(r, theta, gint, ext)
% Total flux function: internal zonal harmonics + external model.
% ext.type = 'ring': ext.a0, ext.sigma, ext.Bc (distributed ring current)
% ext.type = 'harm': ext.G = [G_1..G_Ne] (external zonal harmonics)
Xi = 0;
for n = 1:numel(gint)
    [Phi,~,~] = flux_basis(r, theta, n);
    Xi = Xi + gint(n)*Phi;
end
switch lower(ext.type)
    case 'ring'
        R = r*sin(theta); Z = r*cos(theta);
        Xi = Xi + ring_flux_distributed(R, Z, ext.a0, ext.sigma, ext.Bc);
    case 'harm'
        for n = 1:numel(ext.G)
            [Phi,~,~] = flux_basis_ext(r, theta, n);
            Xi = Xi + ext.G(n)*Phi;
        end
    otherwise
        error('Unknown external type: %s', ext.type);
end
end
