function obs = make_flux_line_mixed(r_start, gint_true, ext, Npts, sigma_r, theta_start)
% Generate a synthetic field line for a mixed (internal + external) flux function.
if nargin < 6, theta_start = pi/2; end
theta_end = asin(sqrt(sin(theta_start)^2 / r_start));
theta = linspace(theta_start, theta_end, Npts)';
fun = @(r,th) flux_mixed(r, th, gint_true, ext);
C = fun(r_start, theta_start);
r_true = zeros(Npts,1);
for i = 1:Npts
    r_guess = r_start*sin(theta(i))^2;
    r_true(i) = trace_level_set_fun(theta(i), C, fun, r_guess, [0.9,3.2]);
end
r_obs = r_true + sigma_r*randn(Npts,1);
r_obs(1) = r_true(1);
[~,~,B_start_true] = field_from_flux(fun, r_start, theta_start);
[~,~,B_end_true]   = field_from_flux(fun, r_true(end), theta(end));
obs.r_start = r_start;
obs.theta = theta;
obs.r_true = r_true;
obs.r_obs = r_obs;
obs.C = C;
obs.B_start_true = B_start_true;
obs.B_end_true = B_end_true;
obs.fun = fun;
end
