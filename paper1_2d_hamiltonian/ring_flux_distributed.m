function Psi = ring_flux_distributed(R, Z, a0, sigma_a, Bc)
% Smooth distributed ring current: Gaussian-weighted sum of circular loops.
% a0, sigma_a in R_E; Bc = total field at the origin (nT).
if nargin < 3, a0 = 4; end
if nargin < 4, sigma_a = 0.8; end
if nargin < 5, Bc = -50; end
a_vals = linspace(max(1.2, a0 - 2*sigma_a), a0 + 2*sigma_a, 9);
w = exp(-0.5*((a_vals - a0)/sigma_a).^2);
w = w / sum(w);
Psi = 0;
for k = 1:numel(a_vals)
    Psi = Psi + w(k) * ring_flux_loop(R, Z, a_vals(k), Bc);
end
end
