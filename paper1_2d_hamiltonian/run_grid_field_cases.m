% Grid-based |B| error for the g1-free cases:
%   ring / no-ring  x  path-only / +endpoint |B|
% 20 realistic 5-year IGRF SV truths, 0.1 km path noise.
% The grid covers the source-free inner-magnetosphere region.
clear; clc;
R_E=6371.2;
g_base=igrf_m0_coeffs(2020);
start_years=1920:5:2015;
nset=numel(start_years);
DG=zeros(10,nset);
for k=1:nset
    g1=igrf_m0_coeffs(start_years(k));
    g2=igrf_m0_coeffs(start_years(k)+5);
    DG(:,k)=g2-g1;
end
Ni=10; Ne=4; sig_km=0.1; lambda=1e-4; Npts=100;
r_starts=[1.2 1.5 1.9 2.4 2.7];
ext_ring=struct('type','ring','a0',6,'sigma',0.8,'Bc',-50);
ext_none=struct('type','harm','G',zeros(4,1));

% fixed spatial grid (r, theta), northern hemisphere
rg=linspace(1.05,3.10,30);
tg=linspace(pi/18,pi/2,30);
[RG,TG]=meshgrid(rg,tg);

cases={'ring_path','ring_B','none_path','none_B'};
ncase=numel(cases);
grid_rms=nan(nset,ncase); grid_p50=nan(nset,ncase);
grid_p95=nan(nset,ncase); grid_max=nan(nset,ncase);
path_B=nan(nset,ncase); path_L=nan(nset,ncase);

for iset=1:nset
    dg_true=DG(:,iset);
    g_truth=g_base+dg_true;

    % observations for the ring truth
    rng(8000+iset); obsR=struct([]);
    for j=1:5
        o=make_flux_line_mixed(r_starts(j),g_truth,ext_ring,Npts,sig_km/R_E);
        o.B_start_obs=o.B_start_true; o.B_end_obs=o.B_end_true;
        obsR=[obsR;o]; %#ok<AGROW>
    end
    % observations for the no-ring truth (same noise seed)
    rng(8000+iset); obsN=struct([]);
    for j=1:5
        o=make_flux_line_mixed(r_starts(j),g_truth,ext_none,Npts,sig_km/R_E);
        o.B_start_obs=o.B_start_true; o.B_end_obs=o.B_end_true;
        obsN=[obsN;o]; %#ok<AGROW>
    end

    % truth grids
    fR=@(r,th) flux_mixed(r,th,g_truth,ext_ring);
    fN=@(r,th) flux_mixed(r,th,g_truth,ext_none);
    BtrueR=grid_bmag(fR,RG,TG);
    BtrueN=grid_bmag(fN,RG,TG);

    % inversions
    [dg_rp,G_rp]=invert_flux_perturbation_full(obsR,g_base,Ni,Ne,lambda);
    [dg_rb,G_rb]=invert_flux_perturbation_B_iter_full(obsR,g_base,Ni,Ne, ...
        lambda,sig_km/R_E,1.0,20);
    [dg_np,G_np]=invert_flux_perturbation_full(obsN,g_base,Ni,Ne,lambda);
    [dg_nb,G_nb]=invert_flux_perturbation_B_iter_full(obsN,g_base,Ni,Ne, ...
        lambda,sig_km/R_E,1.0,20);

    est={dg_rp,G_rp; dg_rb,G_rb; dg_np,G_np; dg_nb,G_nb};
    Btrue={BtrueR,BtrueR,BtrueN,BtrueN};
    obsall={obsR,obsR,obsN,obsN};

    for ic=1:ncase
        dg=est{ic,1}; G=est{ic,2};
        ext_e=struct('type','harm','G',G);
        fe=@(r,th) flux_mixed(r,th,g_base+dg,ext_e);
        Be=grid_bmag(fe,RG,TG);
        err=Be-Btrue{ic};
        grid_rms(iset,ic)=sqrt(mean(err(:).^2));
        grid_p50(iset,ic)=median(abs(err(:)));
        grid_p95(iset,ic)=prctile(abs(err(:)),95);
        grid_max(iset,ic)=max(abs(err(:)));
        [~,~,L]=ensemble_metrics(obsall{ic},dg_true,dg,fe,R_E);
        % along-path |B| RMS recomputed directly
        Bs=[];
        for j=1:numel(obsall{ic})
            o=obsall{ic}(j);
            for i=1:numel(o.r_true)
                [~,~,Bt]=field_from_flux(o.fun,o.r_true(i),o.theta(i));
                [~,~,Be2]=field_from_flux(fe,o.r_true(i),o.theta(i));
                Bs=[Bs;Be2-Bt]; %#ok<AGROW>
            end
        end
        path_B(iset,ic)=sqrt(mean(Bs.^2));
        path_L(iset,ic)=L;
    end
    fprintf('iset %2d done\n',iset);
end

fprintf('\nGrid |B| error (nT), median [16-84] over %d truths:\n',nset);
for ic=1:ncase
    fprintf('%-9s RMS=%7.2f [%6.2f,%7.2f]  p50=%7.2f  p95=%7.2f  max=%8.2f  pathB=%7.2f  pathL=%.3f km\n', ...
        cases{ic},median(grid_rms(:,ic)),prctile(grid_rms(:,ic),16), ...
        prctile(grid_rms(:,ic),84),median(grid_p50(:,ic)), ...
        median(grid_p95(:,ic)),median(grid_max(:,ic)), ...
        median(path_B(:,ic)),median(path_L(:,ic)));
end
save('grid_field_cases.mat','cases','grid_rms','grid_p50','grid_p95','grid_max', ...
    'path_B','path_L','DG','start_years','g_base','rg','tg');

function B=grid_bmag(fun,R,T)
B=zeros(size(R));
for k=1:numel(R)
    [~,~,B(k)]=field_from_flux(fun,R(k),T(k));
end
end
