% SI sensitivity: external truncation Ne and leave-one-out coverage at N=5.
% 20 realistic 5-year IGRF SV truths, free dipole, ring current, +endpoint |B|.
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
Ni=10; lambda=1e-4; sig=0.1/R_E; Npts=100;
r_starts=[1.2 1.5 1.9 2.4 2.7];
ext_ring=struct('type','ring','a0',6,'sigma',0.8,'Bc',-50);
rg=linspace(1.05,3.10,30); tg=linspace(pi/18,pi/2,30);
[RG,TG]=meshgrid(rg,tg);

Ne_list=[2 3 4 5 6]; nNe=numel(Ne_list);
rmse_Ne=nan(nset,nNe); p95_Ne=nan(nset,nNe);
pB_Ne=nan(nset,nNe); pL_Ne=nan(nset,nNe);
rmse_LOO=nan(nset,5); p95_LOO=nan(nset,5);
pB_LOO=nan(nset,5); pL_LOO=nan(nset,5);

for iset=1:nset
    dg_true=DG(:,iset); g_truth=g_base+dg_true;
    rng(8000+iset); obs=struct([]);
    for j=1:5
        o=make_flux_line_mixed(r_starts(j),g_truth,ext_ring,Npts,sig);
        o.B_start_obs=o.B_start_true; o.B_end_obs=o.B_end_true;
        obs=[obs;o]; %#ok<AGROW>
    end
    fT=@(r,th) flux_mixed(r,th,g_truth,ext_ring);
    Btrue=grid_bmag(fT,RG,TG);

    for k=1:nNe
        Ne=Ne_list(k);
        [dg,G]=invert_flux_perturbation_B_iter_full(obs,g_base,Ni,Ne, ...
            lambda,sig,1.0,20);
        [rmse_Ne(iset,k),p95_Ne(iset,k),pB_Ne(iset,k),pL_Ne(iset,k)]= ...
            evalcase(obs,dg_true,dg,G,g_base,Btrue,RG,TG,R_E);
    end
    for j=1:5
        keep=setdiff(1:5,j); ob=obs(keep);
        [dg,G]=invert_flux_perturbation_B_iter_full(ob,g_base,Ni,4, ...
            lambda,sig,1.0,20);
        [rmse_LOO(iset,j),p95_LOO(iset,j),pB_LOO(iset,j),pL_LOO(iset,j)]= ...
            evalcase(ob,dg_true,dg,G,g_base,Btrue,RG,TG,R_E);
    end
    fprintf('truth %2d/%2d done\n',iset,nset);
end

fprintf('\nNe sweep at N=5 (+B):\n');
for k=1:nNe
    fprintf('Ne=%d  gridRMS=%7.2f [%6.2f,%7.2f]  p95=%6.2f  pathB=%6.2f  pathL=%.4f\n', ...
        Ne_list(k),med(rmse_Ne(:,k)),prc(rmse_Ne(:,k),16),prc(rmse_Ne(:,k),84), ...
        med(p95_Ne(:,k)),med(pB_Ne(:,k)),med(pL_Ne(:,k)));
end
fprintf('\nLeave-one-out at N=5 (+B): removed L-shell\n');
for j=1:5
    fprintf('remove %.1f  gridRMS=%7.2f [%6.2f,%7.2f]  p95=%6.2f  pathB=%6.2f  pathL=%.4f\n', ...
        r_starts(j),med(rmse_LOO(:,j)),prc(rmse_LOO(:,j),16),prc(rmse_LOO(:,j),84), ...
        med(p95_LOO(:,j)),med(pB_LOO(:,j)),med(pL_LOO(:,j)));
end
save('sensitivity_structure.mat','Ne_list','rmse_Ne','p95_Ne','pB_Ne','pL_Ne', ...
    'rmse_LOO','p95_LOO','pB_LOO','pL_LOO','r_starts','start_years','g_base','rg','tg');

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
