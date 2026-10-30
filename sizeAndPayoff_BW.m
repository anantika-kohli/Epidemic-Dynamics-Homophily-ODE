%% sizeAndPayoff_BW.m
% Black-and-white version of sizeAndPayoff.m for print.

clear; clc; close all;

OUTFILE = 'Figure_TimeCourse_Payoffs_BW.pdf';

% Figure size in pixels: s = [left bottom width height]. Change this to resize.
s = [100 100 700 320];

%% Panel A: epidemic time course (h = 1) 
PA.beta  = 1.1e-5; PA.c = 8; PA.gamma = 1/5; PA.xi = 0.025; PA.eps = 0.7;
PA.L = 0.7; PA.alpha = 2; PA.kappa = 0.2;
PA.homophily_type = 'power';
PA.pv = 0.05; PA.pa = 0.05;
PA.bV = 0.20; PA.bH = 0.10; PA.bA = 0.15;
PA.cV = 0.08; PA.cH = 0.03; PA.cA = 0.02;
PA.omegaV_min = 0; PA.omegaV_max = 1.5;
PA.omegaH_min = 0; PA.omegaH_max = 0.8;
PA.omegaA_min = 0; PA.omegaA_max = 0.2;
PA.I0 = 200; PA.m = 3;
PA.h = 1;

y0A    = [5000 3999 1000 0 0 1 0 0];   % seed IH = 1
tspanA = [0 500];
[tA, yA] = ode45(@(t,y) homophily_replicator_hill_ode(t,y,PA), tspanA, y0A);
IA_total = yA(:,5) + yA(:,6) + yA(:,7);   % total infected I(t)

%% Panel B: payoffs vs time (h = 1.0) 
PB.beta  = 1.1e-5; PB.c = 8; PB.gamma = 1/5; PB.xi = 0.025; PB.eps = 0.7;
PB.L = 0.7; PB.alpha = 2; PB.kappa = 0.2;
PB.homophily_type = 'assort';
PB.pv = 0.05; PB.pa = 0.05;
PB.bV = 0.20; PB.bH = 0.10; PB.bA = 0.15;
PB.cV = 0.08; PB.cH = 0.03; PB.cA = 0.02;
PB.omegaV_min = 0; PB.omegaV_max = 1.5;
PB.omegaH_min = 0; PB.omegaH_max = 0.8;
PB.omegaA_min = 0; PB.omegaA_max = 0.2;
PB.I0 = 200; PB.m = 3;
PB.h = 1.0;

seedCol = 6;                            % IH
y0B = [5000 4000 1000 0 0 0 0 0];
y0B(seedCol)   = y0B(seedCol) + 1;
y0B(seedCol-4) = y0B(seedCol-4) - 1;
tspanB = [0 500];
tView  = 200;

[tB, yB] = ode45(@(t,y) homophily_replicator_hill_ode(t,y,PB), tspanB, y0B);

n = numel(tB);
PiV = zeros(n,1); PiH = zeros(n,1); PiA = zeros(n,1);
for k = 1:n
    [~, ex] = homophily_replicator_hill_ode(tB(k), yB(k,:).', PB);
    PiV(k) = ex.PiV; PiH(k) = ex.PiH; PiA(k) = ex.PiA;
end

%% Combined figure: (A) and (B) side by side 
fig = figure('Units','pixels','Position',s, 'Color','w');
t = tiledlayout(fig, 1, 2, 'Padding','compact', 'TileSpacing','compact');

% Panel A 
ax1 = nexttile(t);
plot(ax1, tA, IA_total, 'LineWidth', 2, 'Color', 'k');
xlabel(ax1, 'Time'); ylabel(ax1, 'Total infected I(t)');
xlim(ax1, [0 150]); ylim(ax1, [0 600]);
grid(ax1, 'on');
set(ax1, 'Box','off','TickDir','out','LineWidth',0.8,'FontSize',11, ...
    'FontName','Helvetica','XColor','k','YColor','k','GridColor',[0.7 0.7 0.7]);
title(ax1, 'A', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Helvetica');
ax1.TitleHorizontalAlignment = 'left';

% Panel B 
ax2 = nexttile(t);
hold(ax2, 'on');
plot(ax2, tB, PiV, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '-');
plot(ax2, tB, PiH, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '--');
plot(ax2, tB, PiA, 'LineWidth', 2, 'Color', 'k', 'LineStyle', ':');
hold(ax2, 'off');
xlabel(ax2, 'Time'); ylabel(ax2, 'Payoff  \Pi_i');
xlim(ax2, [0 tView]); ylim(ax2, [0 0.15]);
lgd = legend(ax2, {'\Pi_V','\Pi_H','\Pi_A'}, 'Location', 'southeast', ...
    'Box','on', 'FontSize', 9);
title(lgd, 'Payoff');
grid(ax2, 'on');
set(ax2, 'Box','off','TickDir','out','LineWidth',0.8,'FontSize',11, ...
    'FontName','Helvetica','XColor','k','YColor','k','GridColor',[0.7 0.7 0.7]);
title(ax2, 'B', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Helvetica');
ax2.TitleHorizontalAlignment = 'left';

drawnow;

%% Export 
if isAbsolutePath(OUTFILE)
    outPath = OUTFILE;
else
    outPath = fullfile(pwd, OUTFILE);
end

try
    exportgraphics(fig, outPath, 'ContentType', 'vector');
    fprintf('Saved PDF to: %s\n', outPath);
catch ME
    fprintf(2, 'exportgraphics failed: %s\n', ME.message);
end

%% Local function
function tf = isAbsolutePath(pathStr)
tf = ~isempty(regexp(pathStr, '^([A-Za-z]:\\|\\\\|/)', 'once'));
end
