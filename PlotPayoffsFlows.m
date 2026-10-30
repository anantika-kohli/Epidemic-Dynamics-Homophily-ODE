clear; clc; close all;


%% Parameters (edit to match the run you want to inspect) 
P.beta  = 1.0e-5; P.c = 8; P.gamma = 1/5; P.xi = 0.025; P.eps = 0.7;
P.L = 0.7; P.alpha = 2; P.kappa = 0.2;
P.homophily_type = 'assort';        % 'power','assort','logistic'

P.pv = 0.05; P.pa = 0.05;          % switching rates (fixed for this diagnostic)

P.bV = 0.20; P.bH = 0.10; P.bA = 0.15;      % benefits
P.cV = 0.08; P.cH = 0.03; P.cA = 0.02;      % costs
P.omegaV_min = 0; P.omegaV_max = 1.5;
P.omegaH_min = 0; P.omegaH_max = 0.8;
P.omegaA_min = 0; P.omegaA_max = 0.2;
P.I0 = 200; P.m = 3;

%% Initial conditions (seed one infected; default I_H) 
seedCol = 6;                        % 5=IV, 6=IH, 7=IA
y0 = [5000 4000 1000 0 0 0 0 0];
y0(seedCol)   = y0(seedCol) + 1;
y0(seedCol-4) = y0(seedCol-4) - 1;

tspan = [0 500];
tView = 200;                        % x axis limit for viewing

%% Output PDF 
pdfName = 'PayoffsFlows_all.pdf';
if isfile(pdfName); delete(pdfName); end

%% Loop over homophily index 
hvals = 0.1:0.1:1.0;

for hi = 1:numel(hvals)
    P.h = hvals(hi);

    % the model
    [t, y] = ode45(@(t,y) homophily_replicator_hill_ode(t,y,P), tspan, y0);

    % payoffs and flows at each saved time point
    n = numel(t);
    PiV = zeros(n,1); PiH = zeros(n,1); PiA = zeros(n,1);
    THV = zeros(n,1); THA = zeros(n,1);
    for k = 1:n
        [~, ex] = homophily_replicator_hill_ode(t(k), y(k,:).', P);
        PiV(k) = ex.PiV; PiH(k) = ex.PiH; PiA(k) = ex.PiA;
        THV(k) = ex.THV; THA(k) = ex.THA;
    end

    % Figures: one per h 
    figure('Name', sprintf('h = %.1f', P.h), 'Position', [100 100 700 620]);

    % Top: payoffs
    subplot(2,1,1); hold on;
    plot(t, PiV, 'LineWidth', 2.5);
    plot(t, PiH, 'LineWidth', 2.5);
    plot(t, PiA, 'LineWidth', 2.5);
    xlabel('Time'); ylabel('Payoff  \Pi_i');
    title(sprintf('Payoffs vs time   (h = %.1f)', P.h));
    legend({'\Pi_V','\Pi_H','\Pi_A'}, 'Location', 'northeast');
    grid on; xlim([0 tView]); ylim([0 0.15]); styleAxes;

    % Bottom: switching flows
    subplot(2,1,2); hold on;
    plot(t, THV, 'LineWidth', 2.5);
    plot(t, THA, 'LineWidth', 2.5);
    xlabel('Time'); ylabel('Switching flow  T_{H\rightarrow i}');
    title(sprintf('Attitude switching flows vs time   (h = %.1f)', P.h));
    legend({'T_{H\rightarrow V}','T_{H\rightarrow A}'}, 'Location', 'northeast');
    grid on; xlim([0 tView]); ylim([0 12]); styleAxes;

    
    exportgraphics(gcf, pdfName, 'Append', true);
end

fprintf('Saved all %d figures to %s\n', numel(hvals), pdfName);

%% Local function 
function styleAxes()
set(gca,'Box','off','TickDir','out','LineWidth',0.8,'FontSize',12, ...
    'XColor','k','YColor','k');
end