% Can a larger constellation recover from 5 km tracing error?
% N=5/10/15 at sigma_path=5 km, and N=5 with 100/300/1000 points per line.
% Free dipole, ring Bc=-50 nT, +endpoint |B|, 20 IGRF 5-year SV truths.
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
Ni=10; Ne=4; lambda=1e-4; sig=5.0/R_E;
ext_ring=struct('type','ring','a0',6,'sigma',0.8,'Bc',-50);
rg=linspace(1.05,3.10,30); tg=linspace(pi/18,pi/2,30);
[RG,TG]=meshgrid(rg,tg);

r5=[1.2 1.5 1.9 2.4 2.7];
r10=[1.2 1.35 1.5 1.65 1.85 2.05 2.2 2.35 2.5 2.7];
r15=[1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.55 2.7];
rset={r5,r10,r15};
Nlist=[5 10 15]; nN=numel(Nlist);
Nptlist=[100 300 1000]; nP=numel(Nptlist);

rms_N=nan(nset,nN); p95_N=nan(nset,nN); pB_N=nan(nset,nN); pL_N=nan(nset,nN);
rms_P=nan(nset,nP); p95_P=nan(nset,nP); pB_P=nan(nset,nP); pL_P=nan(nset,nP);

for iset=1:nset
    dg_true=DG(:,iset); g_truth=g_base+dg_true;
    fT=@(r,th) flux_mixed(r,th,g_truth,ext_ring);
    Btrue=grid_bmag(fT,RG,TG);

    for k=1:nN
        rs=rset{k}; nl=Nlist(k);
        rng(8000+iset); obs=struct([]);
        for j=1:nl
            o=make_flux_line_mixed(rs(j),g_truth,ext_ring,100,sig);
            o.B_start_obs=o.B_start_true; o.B_end_obs=o.B_end_true;
            obs=[obs;o]; %#ok<AGROW>
        end
        [dg,G]=invert_flux_perturbation_B_iter_full(obs,g_base,Ni,Ne, ...
            lambda,sig,1.0,20);
        [rms_N(iset,k),p95_N(iset,k),pB_N(iset,k),pL_N(iset,k)]= ...
            evalcase(obs,dg_true,dg,G,g_base,Btrue,RG,TG,R_E);
    end

    for k=1:nP
        Npts=Nptlist(k);
        rng(8000+iset); obs=struct([]);
        for j=1:5
            o=make_flux_line_mixed(r5(j),g_truth,ext_ring,Npts,sig);
            o.B_start_obs=o.B_start_true; o.B_end_obs=o.B_end_true;
            obs=[obs;o]; %#ok<AGROW>
        end
        [dg,G]=invert_flux_perturbation_B_iter_full(obs,g_base,Ni,Ne, ...
            lambda,sig,1.0,20);
        [rms_P(iset,k),p95_P(iset,k),pB_P(iset,k),pL_P(iset,k)]= ...
            evalcase(obs,dg_true,dg,G,g_base,Btrue,RG,TG,R_E);
    end
    fprintf('truth %2d/%2d done\n',iset,nset);
end

fprintf('\n5 km path noise, N-line sweep:\n');
for k=1:nN
    fprintf('N=%2d  gridRMS=%8.2f [%7.2f,%8.2f]  p95=%7.2f  pathB=%7.2f  pathL=%.4f\n', ...
        Nlist(k),med(rms_N(:,k)),prc(rms_N(:,k),16),prc(rms_N(:,k),84), ...
        med(p95_N(:,k)),med(pB_N(:,k)),med(pL_N(:,k)));
end
fprintf('\n5 km path noise, N=5, points-per-line sweep:\n');
for k=1:nP
    fprintf('Npts=%4d  gridRMS=%8.2f [%7.2f,%8.2f]  p95=%7.2f  pathB=%7.2f  pathL=%.4f\n', ...
        Nptlist(k),med(rms_P(:,k)),prc(rms_P(:,k),16),prc(rms_P(:,k),84), ...
        med(p95_P(:,k)),med(pB_P(:,k)),med(pL_P(:,k)));
end
save('tracing_accuracy_mitigation.mat','Nlist','Nptlist','rms_N','p95_N','pB_N','pL_N', ...
    'rms_P','p95_P','pB_P','pL_P','r5','r10','r15','start_years','g_base','rg','tg');

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
