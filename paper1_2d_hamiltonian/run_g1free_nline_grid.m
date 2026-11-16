% N=1..5, g1 free, ring current, path-only vs +endpoint |B|.
% Metrics: along-path |B| RMS, field-line RMS, and fixed-grid |B| error.
% 20 realistic 5-year IGRF SV truths, 0.1 km path noise.
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

% fixed grid: 30 x 30 over the source-free inner magnetosphere
rg=linspace(1.05,3.10,30);
tg=linspace(pi/18,pi/2,30);
[RG,TG]=meshgrid(rg,tg);

nN=5; ncase=2; % 1=path, 2=+B
g_rms=nan(nset,nN,ncase); g_p50=nan(nset,nN,ncase);
g_p95=nan(nset,nN,ncase); g_max=nan(nset,nN,ncase);
pB=nan(nset,nN,ncase); pL=nan(nset,nN,ncase);

for iset=1:nset
    dg_true=DG(:,iset);
    g_truth=g_base+dg_true;
    rng(8000+iset); obs=struct([]);
    for j=1:5
        o=make_flux_line_mixed(r_starts(j),g_truth,ext_ring,Npts,sig_km/R_E);
        o.B_start_obs=o.B_start_true; o.B_end_obs=o.B_end_true;
        obs=[obs;o]; %#ok<AGROW>
    end
    fT=@(r,th) flux_mixed(r,th,g_truth,ext_ring);
    Btrue=grid_bmag(fT,RG,TG);

    for nl=1:nN
        [dg,G]=invert_flux_perturbation_full(obs(1:nl),g_base,Ni,Ne,lambda);
        ext_e=struct('type','harm','G',G);
        fe=@(r,th) flux_mixed(r,th,g_base+dg,ext_e);
        [pB(iset,nl,1),pL(iset,nl,1),g_rms(iset,nl,1),g_p50(iset,nl,1), ...
            g_p95(iset,nl,1),g_max(iset,nl,1)]=metrics_grid(obs(1:nl),dg_true,dg,fe,fT,Btrue,RG,TG,R_E);

        try
            [dgb,Gb]=invert_flux_perturbation_B_iter_full(obs(1:nl),g_base,Ni,Ne, ...
                lambda,sig_km/R_E,1.0,20);
            ext_b=struct('type','harm','G',Gb);
            fb=@(r,th) flux_mixed(r,th,g_base+dgb,ext_b);
            [pB(iset,nl,2),pL(iset,nl,2),g_rms(iset,nl,2),g_p50(iset,nl,2), ...
                g_p95(iset,nl,2),g_max(iset,nl,2)]=metrics_grid(obs(1:nl),dg_true,dgb,fb,fT,Btrue,RG,TG,R_E);
        catch
            % leave NaN
        end
    end
    fprintf('iset %2d done\n',iset);
end

fprintf('\nN  path: gridRMS [16-84]  pathB  pathL | +B: gridRMS [16-84]  pathB  pathL\n');
for nl=1:nN
    fprintf('%d  %8.2f [%7.2f,%8.2f] %7.2f %.3f | %8.2f [%7.2f,%8.2f] %7.2f %.3f\n',nl, ...
        med(g_rms(:,nl,1)),prc(g_rms(:,nl,1),16),prc(g_rms(:,nl,1),84), ...
        med(pB(:,nl,1)),med(pL(:,nl,1)), ...
        med(g_rms(:,nl,2)),prc(g_rms(:,nl,2),16),prc(g_rms(:,nl,2),84), ...
        med(pB(:,nl,2)),med(pL(:,nl,2)));
end
save('g1free_nline_grid.mat','g_rms','g_p50','g_p95','g_max','pB','pL', ...
    'DG','start_years','g_base','r_starts','rg','tg');

function [Bpath,Lpath,rms,p50,p95,mx]=metrics_grid(obs,dg_true,dg,fe,fT,Btrue,RG,TG,R_E)
[~,Bpath,Lpath]=ensemble_metrics(obs,dg_true,dg,fe,R_E);
Be=grid_bmag(fe,RG,TG);
err=Be-Btrue;
rms=sqrt(mean(err(:).^2));
p50=median(abs(err(:)));
p95=prctile(abs(err(:)),95);
mx=max(abs(err(:)));
end

function B=grid_bmag(fun,R,T)
B=zeros(size(R));
for k=1:numel(R)
    [~,~,B(k)]=field_from_flux(fun,R(k),T(k));
end
end

function v=med(x), v=median(x,'omitnan'); end
function v=prc(x,p), v=prctile(x(~isnan(x)),p); end
