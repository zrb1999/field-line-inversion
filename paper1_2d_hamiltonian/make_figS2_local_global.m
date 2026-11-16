% Figure S2: local (along-path) versus global (fixed-grid) |B| error.
clear; clc; close all;
C = paper_style();
figdir = fullfile(fileparts(mfilename('fullpath')),'manuscript','figures');
if ~exist(figdir,'dir'), mkdir(figdir); end
S = load('g1free_nline_grid.mat');
N=1:5;
pathB=arrayfun(@(n)median(S.pB(:,n,1)),N);
gridP=arrayfun(@(n)median(S.g_rms(:,n,1)),N);
pathB_B=arrayfun(@(n)median(S.pB(:,n,2)),N);
gridP_B=arrayfun(@(n)median(S.g_rms(:,n,2)),N);

fig=figure('Units','inches','Position',[1 1 5.2 4.6],'Color','w', ...
    'ToolBar','none','MenuBar','none');
ax=axes; hold(ax,'on');
loglog([0.3 400],[0.3 400],'k--','LineWidth',1.0);
h1=loglog(pathB,gridP,'s-','Color',C.gray,'MarkerFaceColor','w','MarkerSize',7,'LineWidth',1.4);
h2=loglog(pathB_B,gridP_B,'o-','Color',C.blue,'MarkerFaceColor','w','MarkerSize',7,'LineWidth',1.8);
for k=1:5
    text(pathB(k)*1.08,gridP(k)*1.08,sprintf('N=%d',N(k)),'Color',C.gray,'FontSize',7.5);
    text(pathB_B(k)*1.08,gridP_B(k)*1.08,sprintf('N=%d',N(k)),'Color',C.blue,'FontSize',7.5);
end
% no-ring N=5 control (30x30 grid)
Gc=load('grid_field_cases.mat');
ic=find(strcmp(Gc.cases,'none_B'));
noRingB=median(Gc.path_B(:,ic)); noRingG=median(Gc.grid_rms(:,ic));
h3=loglog(noRingB,noRingG,'p','Color',C.green,'MarkerFaceColor','w','MarkerSize',10,'LineWidth',1.6);
text(noRingB*1.1,noRingG*0.85,'no ring, N=5','Color',C.green,'FontSize',7.5);
set(ax,'XScale','log','YScale','log');
xlim([0.3 400]); ylim([0.3 400]);
xlabel('along-path |B| RMS (nT)'); ylabel('global grid |B| RMS (nT)');
legend([h1 h2 h3],{'path only','+ endpoint |B|','no-ring control'}, ...
    'Location','northwest','FontSize',7.5,'Box','off');
text(1.1,1.1,'local = global','FontSize',7.5,'Color',C.gray,'Rotation',45);
title('Local versus global field accuracy','FontSize',8.5,'FontWeight','normal');
set(fig,'PaperPositionMode','auto');
exportgraphics(fig,fullfile(figdir,'figS2_local_global.png'),'Resolution',300);
exportgraphics(fig,fullfile(figdir,'figS2_local_global.pdf'),'ContentType','vector');
fprintf('Wrote %s\n',fullfile(figdir,'figS2_local_global.png'));
