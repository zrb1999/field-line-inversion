% Figure 3: ensemble-median spatial |B| error maps for N=1,3,5 (+ endpoint |B|).
% 20 realistic 5-year IGRF SV truths, free dipole, ring current, 30x30 grid.
% A common symmetric-log color scale is used in all panels.
clear; clc; close all;
C = paper_style();
figdir = fullfile(fileparts(mfilename('fullpath')),'manuscript','figures');
if ~exist(figdir,'dir'), mkdir(figdir); end
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
Ni=10; Ne=4; sig=0.1/R_E; lambda=1e-4; Npts=100;
r_starts=[1.2 1.5 1.9 2.4 2.7];
ext_ring=struct('type','ring','a0',6,'sigma',0.8,'Bc',-50);

% fixed grid, identical to the primary metric grid
rg=linspace(1.05,3.10,30); tg=linspace(pi/18,pi/2,30);
[RG,TG]=meshgrid(rg,tg);

Nlist=[1 3 5]; nN=numel(Nlist);
errMaps=nan([size(RG),nset,nN]);
rmsN=nan(nset,nN);

for iset=1:nset
    dg_true=DG(:,iset);
    g_truth=g_base+dg_true;
    rng(8000+iset); obs=struct([]);
    for j=1:5
        o=make_flux_line_mixed(r_starts(j),g_truth,ext_ring,Npts,sig);
        o.B_start_obs=o.B_start_true; o.B_end_obs=o.B_end_true;
        obs=[obs;o]; %#ok<AGROW>
    end
    fT=@(r,th) flux_mixed(r,th,g_truth,ext_ring);
    Btrue=grid_bmag(fT,RG,TG);
    for kk=1:nN
        nl=Nlist(kk);
        [dg,G]=invert_flux_perturbation_B_iter_full(obs(1:nl),g_base,Ni,Ne, ...
            lambda,sig,1.0,20);
        ext_e=struct('type','harm','G',G);
        fe=@(r,th) flux_mixed(r,th,g_base+dg,ext_e);
        err=grid_bmag(fe,RG,TG)-Btrue;
        errMaps(:,:,iset,kk)=err;
        rmsN(iset,kk)=sqrt(mean(err(:).^2));
    end
    fprintf('truth %2d/%2d done\n',iset,nset);
end

% pointwise median signed error and median per-truth RMS
Emed=median(errMaps,3);
rmsMed=median(rmsN,1);

% symmetric-log transform: t = sign(e)*log10(1+|e|/0.5), saturated at 100 nT
T=sign(Emed).*log10(1+abs(Emed)/0.5);
tmax=log10(1+100/0.5);
T=min(max(T,-tmax),tmax);
levels=linspace(-tmax,tmax,35);
tickVals=[-100 -10 -1 0 1 10 100];
tickPos=sign(tickVals).*log10(1+abs(tickVals)/0.5);

fig=figure('Units','inches','Position',[1 1 7.0 2.75],'Color','w');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for kk=1:nN
    nexttile;
    contourf(RG.*sin(TG),RG.*cos(TG),T(:,:,kk),levels,'LineStyle','none');
    hold on; aa=linspace(0,pi/2,100);
    plot(cos(aa),sin(aa),'k','LineWidth',1.6);
    axis equal; xlim([0 3.1]); ylim([0 3.1]);
    xlabel('x (R_E)'); if kk==1, ylabel('z (R_E)'); end
    caxis([-tmax tmax]);
    title(sprintf('(%c) N=%d, median RMS=%.1f nT',char(96+kk),Nlist(kk),rmsMed(kk)),'FontSize',8);
    set(gca,'FontSize',8);
end
colormap(diverging_map(256));
cb=colorbar; cb.Layout.Tile='east';
cb.Ticks=tickPos; cb.TickLabels=compose('%d',tickVals);
cb.Label.String='\Delta|B| (nT)';
set(fig,'PaperPositionMode','auto');
exportgraphics(fig,fullfile(figdir,'fig3_grid_maps.png'),'Resolution',300);
exportgraphics(fig,fullfile(figdir,'fig3_grid_maps.pdf'),'ContentType','vector');
fprintf('Median RMS: N=1 %.2f, N=3 %.2f, N=5 %.2f nT\n',rmsMed);
fprintf('Wrote %s\n',fullfile(figdir,'fig3_grid_maps.png'));
save('fig3_gridmaps_ensemble.mat','errMaps','rmsN','rmsMed','Nlist','start_years', ...
    'g_base','r_starts','rg','tg');

function B=grid_bmag(fun,R,T)
B=zeros(size(R));
for k=1:numel(R)
    [~,~,B(k)]=field_from_flux(fun,R(k),T(k));
end
end
function cmap=diverging_map(n)
if nargin<1, n=256; end
m=floor(n/2);
lower=[linspace(0.10,1,m)',linspace(0.25,1,m)',linspace(0.65,1,m)'];
upper=[linspace(1,0.75,n-m)',linspace(1,0.15,n-m)',linspace(1,0.10,n-m)'];
cmap=[lower;upper];
end
