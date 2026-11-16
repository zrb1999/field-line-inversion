% At 5 km tracing error, does a lower-order internal model reduce noise?
% N=5, ring Bc=-50 nT, +endpoint |B|, 20 truths, Ni=2/4/6/8/10.
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
Ne=4; lambda=1e-4; sig=5.0/R_E; Npts=100;
r_starts=[1.2 1.5 1.9 2.4 2.7];
ext_ring=struct('type','ring','a0',6,'sigma',0.8,'Bc',-50);
rg=linspace(1.05,3.10,30); tg=linspace(pi/18,pi/2,30);
[RG,TG]=meshgrid(rg,tg);
Ni_list=[2 4 6 8 10]; nN=numel(Ni_list);
rms=nan(nset,nN); p95=nan(nset,nN); pB=nan(nset,nN); pL=nan(nset,nN);

for iset=1:nset
    dg_true=DG(:,iset); g_truth=g_base+dg_true;
    fT=@(r,th) flux_mixed(r,th,g_truth,ext_ring);
    Btrue=grid_bmag(fT,RG,TG);
    rng(8000+iset); obs=struct([]);
    for j=1:5
        o=make_flux_line_mixed(r_starts(j),g_truth,ext_ring,Npts,sig);
        o.B_start_obs=o.B_start_true; o.B_end_obs=o.B_end_true;
        obs=[obs;o]; %#ok<AGROW>
    end
    for k=1:nN
        Ni=Ni_list(k); gb=g_base(1:Ni);
        [dg,G]=invert_flux_perturbation_B_iter_full(obs,gb,Ni,Ne, ...
            lambda,sig,1.0,20);
        ext_e=struct('type','harm','G',G);
        gint=gb+dg;
        fe=@(r,th) flux_mixed(r,th,gint,ext_e);
        Be=grid_bmag(fe,RG,TG); err=Be-Btrue;
        rms(iset,k)=sqrt(mean(err(:).^2));
        p95(iset,k)=prctile(abs(err(:)),95);
        [~,pB(iset,k),pL(iset,k)]=ensemble_metrics(obs,dg_true(1:Ni),dg,fe,R_E);
    end
    fprintf('truth %2d/%2d done\n',iset,nset);
end

fprintf('\n5 km tracing, N=5, internal-order sweep:\n');
for k=1:nN
    fprintf('Ni=%2d  gridRMS=%8.2f [%7.2f,%8.2f]  p95=%7.2f  pathB=%7.2f  pathL=%.4f\n', ...
        Ni_list(k),med(rms(:,k)),prc(rms(:,k),16),prc(rms(:,k),84), ...
        med(p95(:,k)),med(pB(:,k)),med(pL(:,k)));
end
save('model_complexity_test.mat','Ni_list','rms','p95','pB','pL','r_starts', ...
    'start_years','g_base','rg','tg');

function B=grid_bmag(fun,R,T)
B=zeros(size(R));
for k=1:numel(R)
    [~,~,B(k)]=field_from_flux(fun,R(k),T(k));
end
end
function v=med(x), v=median(x,'omitnan'); end
function v=prc(x,p), v=prctile(x(~isnan(x)),p); end
