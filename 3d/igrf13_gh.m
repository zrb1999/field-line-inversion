function gh = igrf13_gh(year)
% IGRF-13 Gauss coefficient vector for a decimal year.
% Ordering: g10, g11, h11, g20, g21, h21, g22, h22, g30, ...
% The vector has 195 elements (n = 1..13).
%
% This is the same coefficient table used by MATLAB's igrfmagm.m
% (aeroigrf13dataTable.mat).  For years > 2020 the IGRF-13 secular
% variation column is used to extrapolate.

persistent DATA YEARS
if isempty(DATA)
    S = load(fullfile(matlabroot,'toolbox','aero','aero','aeroigrf13dataTable.mat'));
    DATA = S.data;
    YEARS = 1900:5:2020;
end

if year <= 2020
    idx = find(YEARS <= year, 1, 'last');
    if isempty(idx)
        error('IGRF-13 is valid from 1900 onward (requested %.2f).', year);
    end
    if idx == numel(YEARS)
        gh = DATA(:, idx);
    else
        f = (year - YEARS(idx))/5;
        gh = (1-f)*DATA(:,idx) + f*DATA(:,idx+1);
    end
else
    gh = DATA(:,end-1) + (year - 2020)*DATA(:,end);
end
end
