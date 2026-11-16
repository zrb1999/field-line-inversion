function g = igrf_m0_coeffs(year)
% IGRF-13 zonal (m=0) Gauss coefficients g10..g100 at a decimal year.
% Returns a 10x1 vector in nT.
thisdir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisdir,'..','3d'));
gh = igrf13_gh(year);
g = gh((1:10)'.^2);
end
