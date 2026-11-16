function r = trace_level_set_fun(theta, C, fun, r_guess, r_bounds)
% Solve fun(r,theta)=C on the physical branch.
if nargin < 5 || isempty(r_bounds), r_bounds = [0.9, 3.2]; end
if nargin < 4 || isempty(r_guess), r_guess = 1.1; end
f = @(r) fun(r,theta) - C;
rlo = max(r_bounds(1), r_guess - 0.35);
rhi = min(r_bounds(2), r_guess + 0.35);
flo = f(rlo); fhi = f(rhi);
if flo == 0, r = rlo; return; end
if fhi == 0, r = rhi; return; end
if flo*fhi < 0
    r = fzero(f, [rlo, rhi]);
    return;
end
rg = linspace(r_bounds(1), r_bounds(2), 161);
fg = arrayfun(f, rg);
idx = find(fg(1:end-1).*fg(2:end) <= 0);
if isempty(idx), r = NaN; return; end
rc = 0.5*(rg(idx)+rg(idx+1));
[~,j] = min(abs(rc - r_guess));
r = fzero(f, [rg(idx(j)), rg(idx(j)+1)]);
end
