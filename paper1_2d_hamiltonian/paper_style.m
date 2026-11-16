function C = paper_style()
% Shared publication style for the Paper 1 figures.
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultTextFontName','Arial');
set(groot,'defaultAxesFontSize',8);
set(groot,'defaultTextFontSize',8);
set(groot,'defaultAxesLineWidth',0.9);
set(groot,'defaultLineLineWidth',1.4);
set(groot,'defaultAxesBox','on');
set(groot,'defaultAxesTickDir','out');
set(groot,'defaultAxesTickLength',[0.015 0.015]);
set(groot,'defaultAxesXMinorTick','on');
set(groot,'defaultAxesYMinorTick','on');
C.blue    = [0.0000 0.4470 0.7410];
C.orange  = [0.8500 0.3250 0.0980];
C.green   = [0.4660 0.6740 0.1880];
C.purple  = [0.4940 0.1840 0.5560];
C.yellow  = [0.9290 0.6940 0.1250];
C.cyan    = [0.3010 0.7450 0.9330];
C.gray    = [0.35 0.35 0.35];
C.lightgray = [0.75 0.75 0.75];
end
