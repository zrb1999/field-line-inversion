% Figure 4: N=5 field error versus ring-current intensity, quiet to storm.
% Median and 95th-percentile grid |B| RMS over 20 truths.
clear; clc; close all;
C = paper_style();
figdir = fullfile(fileparts(mfilename('fullpath')),'manuscript','figures');
if ~exist(figdir,'dir'), mkdir(figdir); end

S = load('sensitivity_obs.mat');
E = load('sensitivity_extended.mat');
G = load('grid_field_cases.mat');
ic = find(strcmp(G.cases,'none_B'));
noRing = median(G.grid_rms(:,ic),'omitnan');
noRing95 = median(G.grid_p95(:,ic),'omitnan');

% ring-current cases: no ring plus -10 ... -300 nT
ringBc = [-10 -25 -50 -75 -100 -200 -300];
ringVals = [noRing, median(E.rms_c(:,1),'omitnan'), median(S.rms_c(:,1),'omitnan'), ...
            median(S.rms_c(:,2),'omitnan'), median(S.rms_c(:,3),'omitnan'), ...
            median(E.rms_c(:,2),'omitnan'), median(E.rms_c(:,3),'omitnan'), ...
            median(E.rms_c(:,4),'omitnan')];
ring95 = [noRing95, median(E.p95_c(:,1),'omitnan'), median(S.p95_c(:,1),'omitnan'), ...
          median(S.p95_c(:,2),'omitnan'), median(S.p95_c(:,3),'omitnan'), ...
          median(E.p95_c(:,2),'omitnan'), median(E.p95_c(:,3),'omitnan'), ...
          median(E.p95_c(:,4),'omitnan')];
labels = [{'no ring'}, arrayfun(@(b) sprintf('%d',abs(b)),ringBc, ...
    'UniformOutput',false)];

x = 1:numel(ringVals);
fig = figure('Units','inches','Position',[1 1 5.0 3.5],'Color','w', ...
    'ToolBar','none','MenuBar','none');
ax = axes; hold(ax,'on');
h1 = plot(x,ringVals,'o-','Color',C.blue,'MarkerFaceColor','w', ...
    'MarkerSize',6.5,'LineWidth',1.6);
h2 = plot(x,ring95,'s--','Color',C.orange,'MarkerFaceColor','w', ...
    'MarkerSize',5.8,'LineWidth',1.3);

% reference Bc=-50 nT is the fourth category after "no ring"
xref = 4;
xline(xref,':','Color',C.gray,'LineWidth',0.9);
text(xref+0.08,11.6,'reference','Color',C.gray,'FontSize',7.2);

for k=4:numel(x)
    text(x(k)+0.08,ringVals(k)+0.35,sprintf('%.1f',ringVals(k)), ...
        'Color',C.blue,'FontSize',7.0);
end
text(x(end)-0.1,ring95(end)+0.45,sprintf('%.1f',ring95(end)),'Color',C.orange,'FontSize',7.0, ...
    'HorizontalAlignment','right');

set(ax,'XTick',x,'XTickLabel',labels,'FontSize',8,'XLim',[0.5 numel(x)+0.5], ...
    'YLim',[0 12.5],'Box','on');
xlabel('|B_{c}| (nT): quiet \rightarrow strong storm');
ylabel('grid |B| RMS (nT)');
legend([h1 h2],{'median','95th percentile'},'Location','northwest','Box','off','FontSize',7.5);
title('Ring-current intensity','FontSize',8.5,'FontWeight','normal');

set(fig,'PaperPositionMode','auto');
exportgraphics(fig,fullfile(figdir,'fig4_sensitivity.png'),'Resolution',300);
exportgraphics(fig,fullfile(figdir,'fig4_sensitivity.pdf'),'ContentType','vector');
fprintf('grid: %s\n',mat2str(ringVals,3));
fprintf('p95:  %s\n',mat2str(ring95,3));
fprintf('Wrote %s\n',fullfile(figdir,'fig4_sensitivity.png'));
