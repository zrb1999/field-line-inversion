% Figure 1: concept and method.
% (a) flux-function contours, field lines, anchors, five L-shells
% (b) internal r^-n vs external r^(n+1) radial scaling
% (c) workflow from traced path to coefficients
clear; clc; close all;
C = paper_style();
figdir = fullfile(fileparts(mfilename('fullpath')),'manuscript','figures');
if ~exist(figdir,'dir'), mkdir(figdir); end

R_E = 6371.2;
gint = igrf_m0_coeffs(2020);
r_starts = [1.2 1.5 1.9 2.4 2.7];

% Representative quiet-time IGRF-13 field plus a distributed ring current.
ext = struct('type','ring','a0',6,'sigma',0.8,'Bc',-50);

fig = figure('Units','inches','Position',[1 1 7.0 6.4],'Color','w');

% -------------------------------------------------------------------------
% (a) flux-function contours and field lines
% -------------------------------------------------------------------------
ax1 = subplot(2,2,1); hold(ax1,'on');
rr = linspace(1.0,3.15,100);
tt = linspace(0,pi/2,100);
[R,T] = meshgrid(rr,tt);
X = R.*sin(T); Z = R.*cos(T);
Xi = zeros(size(R));
for k = 1:numel(R)
    Xi(k) = flux_mixed(R(k),T(k),gint,ext);
end
Xi(R < 1.0) = NaN;

% Earth cross-section
aa = linspace(0,pi/2,120);
fill([0 cos(aa) 0],[0 sin(aa) 0],[0.88 0.88 0.88], ...
     'EdgeColor','none','FaceAlpha',1);
plot(cos(aa),sin(aa),'k','LineWidth',2.2);
text(0.48,0.13,'Earth','FontSize',7.5,'Color',[0.25 0.25 0.25]);

% flux-function contours
contour(X,Z,Xi,22,'LineColor',[0.72 0.72 0.72],'LineWidth',0.7);

% five field lines and anchors
cols = [C.blue; C.orange; C.green; C.purple; C.cyan];
hLines = gobjects(numel(r_starts),1);
hAnch  = gobjects(numel(r_starts),1);
hEnd   = gobjects(numel(r_starts),1);
for j = 1:numel(r_starts)
    o = make_flux_line_mixed(r_starts(j),gint,ext,120,0);
    x = o.r_true.*sin(o.theta);
    z = o.r_true.*cos(o.theta);
    hLines(j) = plot(x,z,'-','Color',cols(j,:),'LineWidth',1.7);
    hAnch(j) = plot(r_starts(j),0,'o','MarkerSize',5.5,'LineWidth',1.1, ...
         'MarkerEdgeColor',cols(j,:),'MarkerFaceColor','w');
    hEnd(j) = plot(x(end),z(end),'s','MarkerSize',4.2,'LineWidth',1.0, ...
         'MarkerEdgeColor',cols(j,:),'MarkerFaceColor',cols(j,:));
end
legend([hLines(1) hAnch(1) hEnd(1)],{'field line','anchor','endpoint |B|'}, ...
       'Location','northwest','FontSize',7,'Box','off');
xlabel('x (R_E)'); ylabel('z (R_E)');
xlim([0 3.15]); ylim([0 3.15]); axis equal;
title('(a) Flux-function contours and field lines','FontSize',8.5,'FontWeight','normal');
text(1.55,2.95,sprintf('L = %.1f-%.1f R_E',r_starts(1),r_starts(end)), ...
     'FontSize',7.5,'Color',C.gray);

% -------------------------------------------------------------------------
% (b) radial scaling of internal and external harmonics
% -------------------------------------------------------------------------
ax2 = subplot(2,2,2); hold(ax2,'on');
r = linspace(1.0,3.2,240);
th = pi/3;
cols_n = [C.blue; C.orange];
for kk = 1:2
    n = 2*kk;
    [Pi,~,~] = flux_basis(r,th,n);
    [Pe,~,~] = flux_basis_ext(r,th,n);
    Pi = abs(Pi/Pi(1));
    Pe = abs(Pe/Pe(1));
    loglog(r,Pi,'--','Color',cols_n(kk,:),'LineWidth',1.6);
    loglog(r,Pe,'-','Color',cols_n(kk,:),'LineWidth',1.6);
end
set(gca,'YScale','log');
xlabel('r (R_E)'); ylabel('normalized |\Phi_n|');
xlim([1 3.2]); ylim([5e-3 2e3]);
legend({'n=2 internal','n=2 external','n=4 internal','n=4 external'}, ...
       'Location','northwest','FontSize',7,'Box','off');
text(2.35,0.012,'internal \propto r^{-n}','FontSize',7.5,'Color',C.blue);
text(2.45,80,'external \propto r^{n+1}','FontSize',7.5,'Color',C.orange);
title('(b) Radial scaling','FontSize',8.5,'FontWeight','normal');

% -------------------------------------------------------------------------
% (c) inversion workflow
% -------------------------------------------------------------------------
ax3 = subplot(2,2,[3 4]); hold(ax3,'on');
axis(ax3,'off'); xlim(ax3,[0 1]); ylim(ax3,[0 1]);
labels = {
    sprintf('Traced path\n(r_i, \\theta_i)')
    sprintf('Anchor difference\n\\chi(r_i)=\\chi(r_a)')
    sprintf('Linear system\nA x = b')
    sprintf('Coefficients\n\\delta g_{1-10}, G_{1-4}')
    sprintf('Reconstructed\nB(r, \\theta)')
    };
w = 0.155; gap = 0.035; x0 = 0.012; y0 = 0.36; h = 0.28;
for k = 1:numel(labels)
    x = x0 + (k-1)*(w+gap);
    rectangle('Position',[x y0 w h],'Curvature',0.12, ...
              'EdgeColor',C.blue,'LineWidth',1.1, ...
              'FaceColor',[0.94 0.97 1.00]);
    text(x+w/2,y0+h/2,labels{k},'HorizontalAlignment','center', ...
         'VerticalAlignment','middle','FontSize',7.2,'Interpreter','tex');
    if k < numel(labels)
        x1 = x+w+0.005; x2 = x+w+gap-0.008;
        plot([x1 x2],[y0+h/2 y0+h/2],'k-','LineWidth',1.0);
        plot(x2,y0+h/2,'k>','MarkerSize',4.5,'MarkerFaceColor','k');
    end
end
text(0.5,0.22, ...
    'one linear solve after anchor differencing; no iterative tracing inside the inversion', ...
    'HorizontalAlignment','center','FontSize',7.2,'Color',C.gray);
text(0.012,0.92,'(c) Inversion workflow','FontSize',8.5,'FontWeight','normal');
text(0.5,0.80,'path geometry \rightarrow distributed constraints \rightarrow field coefficients', ...
    'HorizontalAlignment','center','FontSize',7.4,'Color',C.gray);

% -------------------------------------------------------------------------
% export
% -------------------------------------------------------------------------
set(fig,'PaperPositionMode','auto');
exportgraphics(fig,fullfile(figdir,'fig1_concept_method.png'),'Resolution',300);
exportgraphics(fig,fullfile(figdir,'fig1_concept_method.pdf'),'ContentType','vector');
fprintf('Wrote %s\n',fullfile(figdir,'fig1_concept_method.png'));
