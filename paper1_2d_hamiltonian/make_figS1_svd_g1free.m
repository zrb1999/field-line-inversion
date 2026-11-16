% Supplementary figure S1: singular-value spectra for the g1-free,
% ring-current case at N=1,3,5 with endpoint |B|.
clear; clc; close all;
C=paper_style();
figdir=fullfile(fileparts(mfilename('fullpath')),'manuscript','figures');
if ~exist(figdir,'dir'), mkdir(figdir); end
R_E=6371.2;
g_base=igrf_m0_coeffs(2020);
g1=igrf_m0_coeffs(2015); g2=igrf_m0_coeffs(2020);
g_truth=g_base+(g2-g1);
ext=struct('type','ring','a0',6,'sigma',0.8,'Bc',-50);
r_starts=[1.2 1.5 1.9 2.4 2.7]; sig=0.1/R_E; Npts=100;
rng(8020); obs=struct([]);
for j=1:5
    o=make_flux_line_mixed(r_starts(j),g_truth,ext,Npts,sig);
    o.B_start_obs=o.B_start_true; o.B_end_obs=o.B_end_true;
    obs=[obs;o];
end
Nlist=[1 3 5]; cols=[C.gray;C.blue;C.orange];
fig=figure('Units','inches','Position',[1 1 4.4 3.6],'Color','w'); hold on;
for k=1:3
    [~,~,~,info]=invert_flux_perturbation_B_iter_full(obs(1:Nlist(k)),g_base,10,4,1e-4,sig,1.0,20);
    sv=info.svd/info.svd(1);
    semilogy(1:numel(sv),sv,'o-','Color',cols(k,:),'MarkerFaceColor','w','MarkerSize',4.5,'LineWidth',1.4);
    fprintf('N=%d cond=%.1e\n',Nlist(k),info.cond);
end
set(gca,'YScale','log'); xlabel('singular-value index'); ylabel('\sigma_i / \sigma_1');
xlim([0.8 14.2]); ylim([1e-8 1.5]);
legend({'N=1','N=3','N=5'},'Location','southwest','FontSize',7.5,'Box','off');
title('Singular-value spectra (g_1 free, ring current)','FontSize',8.5,'FontWeight','normal');
set(fig,'PaperPositionMode','auto');
exportgraphics(fig,fullfile(figdir,'figS1_svd_g1free.png'),'Resolution',300);
exportgraphics(fig,fullfile(figdir,'figS1_svd_g1free.pdf'),'ContentType','vector');
fprintf('Wrote %s\n',fullfile(figdir,'figS1_svd_g1free.png'));
