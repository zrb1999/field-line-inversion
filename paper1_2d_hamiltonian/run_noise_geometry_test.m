% Does the direction of a 5 km tracing error matter?
% radial (perpendicular), tangential (along-track), and isotropic 2D noise.
% N=5, free dipole, ring Bc=-50 nT, +endpoint |B|, 20 truths.
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
Ni=10; Ne=4; lambda=1e-4; sig=5.0/R_E; Npts=100;
r_starts=[1.2 1.5 1.9 2.4 2.7];
ext_ring=struct('type','ring','a0',6,'sigma',0.8,'Bc',-50);
rg=linspace(1.05,3.10,30); tg=linspace(pi/18,pi/2,30);
[RG,TG]=meshgrid(rg,tg);
cases={'radial','tangential','isotropic','correlated'}; nc=numel(cases);
rms=nan(nset,nc); p95=nan(nset,nc); pB=nan(nset,nc); pL=nan(nset,nc);

for iset=1:nset
    dg_true=DG(:,iset); g_truth=g_base+dg_true;
    fT=@(r,th) flux_mixed(r,th,g_truth,ext_ring);
    Btrue=grid_bmag(fT,RG,TG);
    rng(8000+iset);
    for ic=1:nc
        obs=struct([]);
        for j=1:5
            o=make_flux_line_mixed(r_starts(j),g_truth,ext_ring,Npts,0);
            th=o.theta; r=o.r_true;
            x=r.*sin(th); z=r.*cos(th);
            if strcmp(cases{ic},'radial')
                rn=r+sig*randn(Npts,1); thn=th;
            elseif strcmp(cases{ic},'tangential')
                tx=gradient(x); tz=gradient(z);
                tn=hypot(tx,tz); tx=tx./tn; tz=tz./tn;
                e=sig*randn(Npts,1);
                xn=x+e.*tx; zn=z+e.*tz;
                rn=hypot(xn,zn); thn=atan2(xn,zn);
            elseif strcmp(cases{ic},'correlated')
                e0=randn(Npts,1);
                k=exp(-0.5*(((1:21)-11)/5).^2); k=k/sum(k);
                e=conv(e0,k,'same'); e=sig*e/std(e);
                rn=r+e; thn=th;
            else
                xn=x+sig*randn(Npts,1); zn=z+sig*randn(Npts,1);
                rn=hypot(xn,zn); thn=atan2(xn,zn);
            end
            rn(1)=r(1); thn(1)=th(1);
            o.r_obs=rn; o.theta=thn;
            o.B_start_obs=o.B_start_true; o.B_end_obs=o.B_end_true;
            obs=[obs;o]; %#ok<AGROW>
        end
        [dg,G]=invert_flux_perturbation_B_iter_full(obs,g_base,Ni,Ne, ...
            lambda,sig,1.0,20);
        [rms(iset,ic),p95(iset,ic),pB(iset,ic),pL(iset,ic)]= ...
            evalcase(obs,dg_true,dg,G,g_base,Btrue,RG,TG,R_E);
    end
    fprintf('truth %2d/%2d done\n',iset,nset);
end

fprintf('\n5 km tracing error, noise direction (N=5, +B):\n');
for ic=1:nc
    fprintf('%-11s gridRMS=%8.2f [%7.2f,%8.2f]  p95=%7.2f  pathB=%7.2f  pathL=%.4f\n', ...
        cases{ic},med(rms(:,ic)),prc(rms(:,ic),16),prc(rms(:,ic),84), ...
        med(p95(:,ic)),med(pB(:,ic)),med(pL(:,ic)));
end
save('noise_geometry_test.mat','cases','rms','p95','pB','pL','r_starts', ...
    'start_years','g_base','rg','tg');

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
