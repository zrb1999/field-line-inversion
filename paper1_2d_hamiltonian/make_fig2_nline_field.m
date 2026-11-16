% Figure 2: N-line field metrics with free dipole and ring current.
% (a) global grid |B| RMS, (b) along-path |B| RMS,
% (c) field-line RMS, (d) grid 95th percentile.
clear; clc; close all;
C = paper_style();
figdir = fullfile(fileparts(mfilename('fullpath')),'manuscript','figures');
if ~exist(figdir,'dir'), mkdir(figdir); end
S = load('g1free_nline_grid.mat');
N = (1:5)';
med = @(x) median(x,'omitnan');
lo  = @(x) prctile(x(~isnan(x)),16);
hi  = @(x) prctile(x(~isnan(x)),84);

fig = figure('Units','inches','Position',[1 1 7.0 5.4],'Color','w');

% (a) global grid RMS
ax1=subplot(2,2,1); hold(ax1,'on');
plot_band(N,S.g_rms(:,:,1),C.gray);
plot_band(N,S.g_rms(:,:,2),C.blue);
semilogy(N,arrayfun(@(n)med(S.g_rms(:,n,1)),N),'o-','Color',C.gray,'MarkerFaceColor','w','MarkerSize',5.5);
semilogy(N,arrayfun(@(n)med(S.g_rms(:,n,2)),N),'o-','Color',C.blue,'MarkerFaceColor','w','MarkerSize',5.5);
set(ax1,'YScale','log'); xlim([0.8 5.2]); xticks(1:5); ylim([1 400]);
ylabel('global grid |B| RMS (nT)'); xlabel('number of field lines');
text(3.55,2.2,'+B','Color',C.blue,'FontSize',8,'FontWeight','bold');
text(3.55,250,'path only','Color',C.gray,'FontSize',8,'FontWeight','bold');
title('(a) Global field error','FontSize',8.5,'FontWeight','normal');

% (b) along-path |B| RMS
ax2=subplot(2,2,2); hold(ax2,'on');
plot_band(N,S.pB(:,:,1),C.gray);
plot_band(N,S.pB(:,:,2),C.blue);
semilogy(N,arrayfun(@(n)med(S.pB(:,n,1)),N),'s-','Color',C.gray,'MarkerFaceColor','w','MarkerSize',5.5);
semilogy(N,arrayfun(@(n)med(S.pB(:,n,2)),N),'s-','Color',C.blue,'MarkerFaceColor','w','MarkerSize',5.5);
set(ax2,'YScale','log'); xlim([0.8 5.2]); xticks(1:5); ylim([0.3 500]);
ylabel('along-path |B| RMS (nT)'); xlabel('number of field lines');
text(3.55,0.45,'+B','Color',C.blue,'FontSize',8,'FontWeight','bold');
text(3.55,250,'path only','Color',C.gray,'FontSize',8,'FontWeight','bold');
title('(b) Along-path field error','FontSize',8.5,'FontWeight','normal');

% (c) field-line RMS
ax3=subplot(2,2,3); hold(ax3,'on');
plot_band(N,S.pL(:,:,1),C.gray);
plot_band(N,S.pL(:,:,2),C.blue);
plot(N,arrayfun(@(n)med(S.pL(:,n,1)),N),'^-','Color',C.gray,'MarkerFaceColor','w','MarkerSize',5.5);
plot(N,arrayfun(@(n)med(S.pL(:,n,2)),N),'^-','Color',C.blue,'MarkerFaceColor','w','MarkerSize',5.5);
xlim([0.8 5.2]); xticks(1:5); ylim([0 0.11]);
ylabel('field-line RMS (km)'); xlabel('number of field lines');
text(3.65,0.100,'path only','Color',C.gray,'FontSize',8,'FontWeight','bold');
text(3.65,0.034,'+B','Color',C.blue,'FontSize',8,'FontWeight','bold');
title('(c) Path reconstruction','FontSize',8.5,'FontWeight','normal');

% (d) grid p95
ax4=subplot(2,2,4); hold(ax4,'on');
plot_band(N,S.g_p95(:,:,1),C.gray);
plot_band(N,S.g_p95(:,:,2),C.blue);
semilogy(N,arrayfun(@(n)med(S.g_p95(:,n,1)),N),'d-','Color',C.gray,'MarkerFaceColor','w','MarkerSize',5.5);
semilogy(N,arrayfun(@(n)med(S.g_p95(:,n,2)),N),'d-','Color',C.blue,'MarkerFaceColor','w','MarkerSize',5.5);
set(ax4,'YScale','log'); xlim([0.8 5.2]); xticks(1:5); ylim([1 1000]);
ylabel('global grid |B| p95 (nT)'); xlabel('number of field lines');
text(3.55,2.0,'+B','Color',C.blue,'FontSize',8,'FontWeight','bold');
text(3.55,500,'path only','Color',C.gray,'FontSize',8,'FontWeight','bold');
title('(d) Off-path tail','FontSize',8.5,'FontWeight','normal');

set(fig,'PaperPositionMode','auto');
axs = findall(fig,'Type','axes');
for k = 1:numel(axs)
    if isprop(axs(k),'Toolbar') && ~isempty(axs(k).Toolbar)
        axs(k).Toolbar.Visible = 'off';
    end
end
exportgraphics(fig,fullfile(figdir,'fig2_nline_field.png'),'Resolution',300);
exportgraphics(fig,fullfile(figdir,'fig2_nline_field.pdf'),'ContentType','vector');
fprintf('Wrote %s\n',fullfile(figdir,'fig2_nline_field.png'));

function plot_band(N,X,col)
lo=arrayfun(@(n)prctile(X(~isnan(X(:,n)),n),16),N);
hi=arrayfun(@(n)prctile(X(~isnan(X(:,n)),n),84),N);
fill([N; flipud(N)],[lo; flipud(hi)],col,'FaceAlpha',0.15,'EdgeColor','none');
end
