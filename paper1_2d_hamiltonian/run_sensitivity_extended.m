% Extended sensitivity: large path noise (2, 5 km) and ring current from
% quiet to strong storm (Bc = -10, -100, -200, -300 nT).
% N=5, free dipole, +endpoint |B|, 20 IGRF 5-year SV truths.
clear; clc;
R_E=6371.2;
g_base=igrf_m0_coeffs(2020);
start_years=1920:5:2015; nset=numel(start_years);
DG=zeros(10,nset);
for k=1:nset
    g1=igrf_m0_coeffs(start_years(k));
    g2=igrf_m0_coeffs(start_years(k)+5);
    DG(:,k)=g2-g1;
end
Ni=10; Ne=4; lambda=1e-4; Npts=100;
r_starts=[1.2 1.5 1.9 2.4 2.7];
rg=linspace(1.05,3.10,30); tg=linspace(pi/18,pi/2,30);
[RG,TG]=meshgrid(rg,tg);

sig_list=[2.0 5.0]; ns=numel(sig_list);
Bc_list=[-10 -100 -200 -300]; nc=numel(Bc_list);
rms_s=nan(nset,ns); p95_s=nan(nset,ns); pB_s=nan(nset,ns); pL_s=nan(nset,ns);
rms_c=nan(nset,nc); p95_c=nan(nset,nc); pB_c=nan(nset,nc); pL_c=nan(nset,nc);

for iset=1:nset
    dg_true=DG(:,iset); g_truth=g_base+dg_true;

    % large path noise
    for k=1:ns
        sig=sig_list(k)/R_E;
        ext_ring=struct('type','ring','a0',6,'sigma',0.8,'Bc',-50);
        fT=@(r,th) flux_mixed(r,th,g_truth,ext_ring);
        Btrue=grid_bmag(fT,RG,TG);
        rng(8000+iset); obs=struct([]);
        for j=1:5
            o=make_flux_line_mixed(r_starts(j),g_truth,ext_ring,Npts,sig);
            o.B_start_obs=o.B_start_true; o.B_end_obs=o.B_end_true;
            obs=[obs;o]; %#ok<AGROW>
        end
        [dg,G]=invert_flux_perturbation_B_iter_full(obs,g_base,Ni,Ne, ...
            lambda,sig,1.0,20);
        [rms_s(iset,k),p95_s(iset,k),pB_s(iset,k),pL_s(iset,k)]= ...
            evalcase(obs,dg_true,dg,G,g_base,Btrue,RG,TG,R_E);
    end

    % ring current from quiet to strong storm
    for k=1:nc
        ext_c=struct('type','ring','a0',6,'sigma',0.8,'Bc',Bc_list(k));
        fTc=@(r,th) flux_mixed(r,th,g_truth,ext_c);
        BtrueC=grid_bmag(fTc,RG,TG);
        sig=0.1/R_E;
        rng(8000+iset); obs=struct([]);
        for j=1:5
            o=make_flux_line_mixed(r_starts(j),g_truth,ext_c,Npts,sig);
            o.B_start_obs=o.B_start_true; o.B_end_obs=o.B_end_true;
            obs=[obs;o]; %#ok<AGROW>
        end
        [dg,G]=invert_flux_perturbation_B_iter_full(obs,g_base,Ni,Ne, ...
            lambda,sig,1.0,20);
        [rms_c(iset,k),p95_c(iset,k),pB_c(iset,k),pL_c(iset,k)]= ...
            evalcase(obs,dg_true,dg,G,g_base,BtrueC,RG,TG,R_E);
    end
    fprintf('truth %2d/%2d done\n',iset,nset);
end

fprintf('\nLarge path-noise cases (N=5, +B, ring Bc=-50 nT):\n');
for k=1:ns
    fprintf('sigma=%.1f km  gridRMS=%8.2f [%7.2f,%8.2f]  p95=%7.2f  pathB=%7.2f  pathL=%.4f\n', ...
        sig_list(k),med(rms_s(:,k)),prc(rms_s(:,k),16),prc(rms_s(:,k),84), ...
        med(p95_s(:,k)),med(pB_s(:,k)),med(pL_s(:,k)));
end
fprintf('\nRing-current sweep (N=5, +B, 0.1 km):\n');
for k=1:nc
    fprintf('Bc=%5.0f nT  gridRMS=%8.2f [%7.2f,%8.2f]  p95=%7.2f  pathB=%7.2f  pathL=%.4f\n', ...
        Bc_list(k),med(rms_c(:,k)),prc(rms_c(:,k),16),prc(rms_c(:,k),84), ...
        med(p95_c(:,k)),med(pB_c(:,k)),med(pL_c(:,k)));
end
save('sensitivity_extended.mat','sig_list','Bc_list','rms_s','p95_s','pB_s','pL_s', ...
    'rms_c','p95_c','pB_c','pL_c','r_starts','start_years','g_base','rg','tg');

function [rms,p95,pB,pL]=evalcase(obs,dg_true,dg,G,g_base,Btrue,RG,TG,R_E)
ext_e=struct('type','harm','G',G);
gint=g_base+dg;
fe=@(r,th) flux_mixed(r,th,gint,ext_e);
Be=grid_bmag(fe,RG,TG); err=Be-Btrue;
rms=sqrt(mean(err(:).^2)); p95=prctile(abs(err(:)),95);
[~,pB,pL]=ensemble_metrics(obs,dg_true,dg,fe,R_E);
end

function B=grid_bmag(fun,R,T)
B=zeros(size(R));
for k=1:numel(R)
    [~,~,B(k)]=field_from_flux(fun,R(k),T(k));
end
end
function v=med(x), v=median(x,'omitnan'); end
function v=prc(x,p), v=prctile(x(~isnan(x)),p); end
