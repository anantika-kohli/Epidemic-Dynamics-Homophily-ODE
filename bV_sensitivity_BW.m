%% bV_sensitivity_BW.m
% Black-and-white version of bV_sensitivity.m for print.


clear; clc; close all;


OUTFILE = 'Figure_Sensitivity_bV_BW.pdf';

% Figure size in pixels: s = [left bottom width height]. Change this to resize.
s = [100 100 700 320];

%% Shared parameters 
Pbase.beta  = 1.0e-5; Pbase.c = 8; Pbase.gamma = 1/5; Pbase.xi = 0.025; Pbase.eps = 0.7;
Pbase.L = 0.7; Pbase.alpha = 2; Pbase.kappa = 0.2;
Pbase.cV = 0.03; Pbase.cH = 0.03; Pbase.cA = 0.03;
Pbase.omegaV_min = 0; Pbase.omegaV_max = 1.5;
Pbase.omegaH_min = 0; Pbase.omegaH_max = 0.8;
Pbase.omegaA_min = 0; Pbase.omegaA_max = 0.2;
Pbase.I0 = 200; Pbase.m = 3;

N0    = 10000;
tspan = [0 500];
hvals = linspace(0,1,80);
seedCols   = [5 6 7];
seedNames  = {'pro-vax','hesitant','anti-vax'};
lineStyles = {'-','--',':'};                       % solid, dashed, dotted

pv = 0.30; pa = 0.30;   % strong adaptation, both panels

%% Panel A: baseline (bV = bH = bA = 0.10)
PA = Pbase;
PA.bV = 0.10; PA.bH = 0.10; PA.bA = 0.10;
covA = run_scenario(PA, 'power', pv, pa, hvals, seedCols, tspan, N0);

%% Panel B: increased pro-vax benefit (bV = 0.15) 
PB = Pbase;
PB.bV = 0.15; PB.bH = 0.10; PB.bA = 0.10;
covB = run_scenario(PB, 'power', pv, pa, hvals, seedCols, tspan, N0);

%% Combined figure: (A) and (B) side by side 
fig = figure('Units','pixels','Position',s, 'Color','w');
t = tiledlayout(fig, 1, 2, 'Padding','compact', 'TileSpacing','compact');

panelData   = {covA, covB};
panelLetter = {'A', 'B'};

for p = 1:2
    ax = nexttile(t);
    hold(ax, 'on');
    for k = 1:3
        plot(ax, hvals, panelData{p}(:,k), 'LineWidth', 2, 'Color', 'k', ...
            'LineStyle', lineStyles{k});
    end
    hold(ax, 'off');

    ylim(ax, [20 90]);
    xlim(ax, [0 1]);
    xlabel(ax, 'Homophily index h');
    if p == 1
        ylabel(ax, 'Final Vaccination Coverage (%)');
    end

    lgd = legend(ax, seedNames, 'Location', 'best', 'Box', 'on', 'FontSize', 9);
    title(lgd, 'First infection seed');

    grid(ax, 'on');
    set(ax, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.8, 'FontSize', 11, ...
        'FontName', 'Helvetica', 'XColor', 'k', 'YColor', 'k', 'XTick', 0:0.2:1, ...
        'GridColor', [0.7 0.7 0.7]);

    title(ax, panelLetter{p}, 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Helvetica');
ax.TitleHorizontalAlignment = 'left';
end

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
