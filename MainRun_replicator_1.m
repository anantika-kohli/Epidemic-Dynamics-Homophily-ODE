clear; clc; close all;
% Driver for the full Section 1.1 replicator model with prevalence-dependent
% omega_i(I). Keep this file and homophily_replicator_hill_ode.m in the same
% folder and run this script.

%% ---------------- Epidemiological parameters ----------------
P.beta  = 1.1e-5;     % transmission probability per contact
P.c     = 8;        % contacts per unit time
P.gamma = 1/5;      % recovery rate (5-day infectious period)
P.xi    = 0.025;    % vaccination rate of pro-vax susceptibles
P.eps   = 0.7;      % vaccine efficacy (breakthrough fraction = 1-eps)

%% ---------------- Homophily ----------------
P.L     = 0.7;      % baseline within-group contact floor
P.alpha = 2;        % power-law exponent
P.kappa = 0.2;      % assortativity curvature (for 'assort')
P.homophily_type = 'power';   % 'power' | 'assort' | 'logistic'
% For 'logistic' also set:  P.k = 10;  P.h0 = 0.5;

%% ---------------- Replicator switching rates (Section 1.1.3) ----------------
% Per Prof. Moghadas: two direction-specific rates instead of a single rho.
P.pv = 0.05;        % rate of hesitant -> pro-vax switching
P.pa = 0.05;        % rate of hesitant -> anti-vax switching (set ~= pv for asymmetry)
P.bV = 0.20; P.bH = 0.10; P.bA = 0.15;   % baseline benefits
P.cV = 0.08; P.cH = 0.03; P.cA = 0.02;   % perceived costs

%% ---------------- Prevalence-dependent omega_i(I) (Section 1.1.4) ----------------
% The document gives the Hill FORM but not the values. These defaults let the
% sensitivity ramp from omega_min (low prevalence) up to omega_max (high
% prevalence). Tune freely. To make the feedback bite harder, LOWER I0,
% RAISE m, or RAISE the omega_max values.
P.omegaV_min = 0;  P.omegaV_max = 1.5;
P.omegaH_min = 0;  P.omegaH_max = 0.8;
P.omegaA_min = 0;  P.omegaA_max = 0.2;
P.I0 = 200;         % half-saturation prevalence (in head-counts)
P.m  = 3;           % Hill steepness

%% ---------------- Initial conditions ----------------
% 50% pro-vax, 40% hesitant, 10% anti-vax; one hesitant individual infected
N0    = 10000;
y0    = [5000 3999 1000 0 0 1 0 0];   % [SV SH SA V IV IH IA R]
tspan = [0 500];

%% ---------------- Sweep over the homophily index ----------------
hvals          = linspace(0,1,100);
attack         = zeros(size(hvals));
cov_stock      = zeros(size(hvals));
cov_cumulative = zeros(size(hvals));

for k = 1:numel(hvals)
    P.h = hvals(k);
    [t,y] = ode45(@(t,y) homophily_replicator_hill_ode(t,y,P), tspan, y0);

    attack(k)         = 100*y(end,8)/N0;                 % R(end): cumulative infected
    cov_stock(k)      = 100*y(end,4)/N0;                 % V(end): currently-protected
    cov_cumulative(k) = 100*trapz(t, P.xi*y(:,1))/N0;    % integral of xi*SV: ever-vaccinated
end

%% ---------------- Figure 1: epidemic size ----------------
figure;
plot(hvals, attack, 'LineWidth', 3);
xlabel('Homophily index h'); ylabel('Attack Rate (%)');
title('Effect of homophily on epidemic size'); ylim([25 50]); grid on; styleAxes;

%% ---------------- Figure 2: vaccination coverage (both metrics) ----------------
figure; hold on;
plot(hvals, cov_cumulative, 'LineWidth', 3, 'DisplayName', 'Total ever vaccinated');
plot(hvals, cov_stock,      'LineWidth', 3, 'DisplayName', 'V at final time (stock)');
xlabel('Homophily index h'); ylabel('Vaccination Coverage (%)'); ylim([40 70]);
title('Vaccination coverage'); legend('Location','best'); grid on; styleAxes;

%% ---------------- Figure 3: epidemic time course at one h ----------------
P.h = 1;
[t,y] = ode45(@(t,y) homophily_replicator_hill_ode(t,y,P), tspan, y0);
figure;
plot(t, y(:,5)+y(:,6)+y(:,7), 'LineWidth', 3);
xlabel('Time'); ylabel('Total infected I(t)');
title('Epidemic time course (h = 1)'); grid on;
set(gca,'Box','off','TickDir','out','LineWidth',0.8,'FontSize',13,'XColor','k','YColor','k');
ylim([0 600]); xlim([0 150]);

%% ---------------- Local function ----------------
function styleAxes()
set(gca,'Box','off','TickDir','out','LineWidth',0.8,'FontSize',13, ...
    'XColor','k','YColor','k','XTick',0:0.1:1);
xlim([0 1]);
end
