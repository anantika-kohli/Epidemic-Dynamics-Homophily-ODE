%% make_panels_BW.m
% Black-and-white version of make_panels.m for print.

clear; clc; close all;

%% Base parameters 
P.beta  = 1.0e-5; P.c = 8; P.gamma = 1/5; P.xi = 0.025; P.eps = 0.7;
P.L = 0.7; P.alpha = 2; P.kappa = 0.2;
P.bV = 0.20; P.bH = 0.05; P.bA = 0.15;
P.cV = 0.08; P.cH = 0.03; P.cA = 0.02;
P.omegaV_min = 0; P.omegaV_max = 1.5;
P.omegaH_min = 0; P.omegaH_max = 0.8;
P.omegaA_min = 0; P.omegaA_max = 0.2;
P.I0 = 200; P.m = 3;

N0    = 10000;
tspan = [0 500];
hvals = linspace(0,1,80);

seedCols  = [5 6 7];                            % columns for IV, IH, IA
seedNames = {'pro-vax','hesitant','anti-vax'};
lineStyles = {'-','--',':'};                    % solid, dashed, dotted

%% FIGURES

PANELS = {
    'A', 'power', 0.30,  0.30, 'coverage'
    'B', 'assort', 0.30,  0.30, 'coverage'
    'C', 'power', 0.30,  0.30, 'infections'
    'D', 'assort', 0.30,  0.30, 'infections'
};

OUTFILE = 'C:\Users\Anantika\OneDrive - York University\Documents\Desktop\ENURA\Figure4ABCD.pdf';

YLIM_COVERAGE   = [20 90];
YLIM_INFECTIONS = [0 6000];


nPanels = size(PANELS,1);

if nPanels <= 2
    nRows = 1; nCols = nPanels;
else
    nCols = 2; nRows = ceil(nPanels/2);
end

fig = figure('Units','inches','Position',[1 1 4.5*nCols 4*nRows], 'Color','w');
t = tiledlayout(fig, nRows, nCols, 'Padding','compact', 'TileSpacing','compact');

for p = 1:nPanels
    label          = PANELS{p,1};
    homophily_type = PANELS{p,2};
    pv             = PANELS{p,3};
    pa             = PANELS{p,4};
    metric         = PANELS{p,5};

    [cov, siz] = run_scenario(P, homophily_type, pv, pa, hvals, seedCols, tspan, N0);

    if strcmp(metric, 'coverage')
        data = cov;
        ylab = 'Final Vaccination Coverage (%)';
        ylims = YLIM_COVERAGE;
    else
        data = siz;
        ylab = 'Total Infections';
        ylims = YLIM_INFECTIONS;
    end

    ax = nexttile(t);
    hold(ax, 'on');
    for s = 1:3
        plot(ax, hvals, data(:,s), 'LineWidth', 2, 'Color', 'k', ...
            'LineStyle', lineStyles{s});
    end
    hold(ax, 'off');

    if ~isempty(ylims)
        ylim(ax, ylims);
    end
    xlim(ax, [0 1]);
    xlabel(ax, 'Homophily index h');
    if mod(p-1, nCols) == 0                 % only label y-axis on leftmost column
        ylabel(ax, ylab);
    end

    % Legend 
    lgd = legend(ax, seedNames, 'Location', 'best', 'Box', 'on', 'FontSize', 9);
    title(lgd, 'First infection seed');

    grid(ax, 'on');
    set(ax, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.8, 'FontSize', 11, ...
        'FontName', 'Helvetica', 'XColor', 'k', 'YColor', 'k', 'XTick', 0:0.2:1, ...
        'GridColor', [0.7 0.7 0.7]);

      % Panel letter 
    title(ax, label, 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Helvetica');
    ax.TitleHorizontalAlignment = 'left';
end

drawnow;

%% Size and export 
% Adjust s = [left bottom width height] in pixels to change the figure size.
set(gcf, 'Units', 'pixels');   
s = [100 100 700 600];
set(gcf, 'Position', s)

exportgraphics(gcf, OUTFILE, 'ContentType', 'vector');
fprintf('Saved PDF to: %s\n', OUTFILE);
