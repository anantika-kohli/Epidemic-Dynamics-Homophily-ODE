function dydt = homophily_replicator_ode(~, y, p)

% State variables
SV = y(1);
SH = y(2);
SA = y(3);
V  = y(4);
IV = y(5);
IH = y(6);
IA = y(7);
R  = y(8);

%% Homophily function
switch lower(p.homophily_type)
    case 'power'
        F = p.L + (1 - p.L) * p.h^p.alpha;

    case 'assort'
        F = p.L + (1 - p.L) * (p.h / (p.h + p.kappa*(1 - p.h) + eps));

    otherwise
        error('Unknown homophily_type. Use ''power'' or ''assort''.');
end

%% Forces of infection
lambdaV = p.beta * p.c * ( F*IV + (1-F)*(IH + IA) );
lambdaH = p.beta * p.c * ( F*IH + (1-F)*(IV + IA) );
lambdaA = p.beta * p.c * ( F*IA + (1-F)*(IV + IH) );

%% Attitude frequencies (excluding recovered)
M = SV + SH + SA + V + IV + IH + IA;

if M > 0
    xV = (SV + V + IV) / M;
    xH = (SH + IH)      / M;
    xA = (SA + IA)      / M;
else
    xV = 0; xH = 0; xA = 0;
end

%% Payoffs
% Pi_i = b_i - c_i - omega_i * lambda_i
PiV = p.bV - p.cV - p.omegaV * lambdaV;
PiH = p.bH - p.cH - p.omegaH * lambdaH;
PiA = p.bA - p.cA - p.omegaA * lambdaA;

%% Replicator-inspired switching terms
% Only hesitant individuals change attitude
THV = p.rho * p.c * (1-F) * SH * xV * max(PiV - PiH, 0);
THA = p.rho * p.c * (1-F) * SH * xA * max(PiA - PiH, 0);

%% ODE system
dSV = -lambdaV * SV + THV - p.xi * SV;
dSH = -lambdaH * SH - THV - THA;
dSA = -lambdaA * SA + THA;

dV  = p.xi * SV - (1 - p.eps) * lambdaV * V;

dIV = lambdaV * SV + (1 - p.eps) * lambdaV * V - p.gamma * IV;
dIH = lambdaH * SH - p.gamma * IH;
dIA = lambdaA * SA - p.gamma * IA;

dR  = p.gamma * (IV + IH + IA);

dydt = [dSV; dSH; dSA; dV; dIV; dIH; dIA; dR];
end
clear; clc;

%% Parameters
params.beta  = 0.00001;
params.c     = 8;
params.gamma = 1/5;
params.xi    = 0.025;
params.eps   = 0.7;

% Homophily parameters
params.L     = 0.7;
params.alpha = 2;      % power-law homophily exponent
params.kappa = 0.2;    % for assortativity-based option if used later

% Replicator / attitude adaptation
params.rho   = 0.02;   % speed of attitude adaptation

% Payoff structure: Pi_i = b_i - c_i - omega_i * lambda_i
params.bV = 0.20;
params.bH = 0.10;
params.bA = 0.15;

params.cV = 0.08;      % perceived cost of pro-vax attitude / vaccination
params.cH = 0.03;      % perceived cost of remaining hesitant
params.cA = 0.02;      % perceived cost of anti-vax attitude

params.omegaV = 1.5;   % weight on infection risk
params.omegaH = 0.8;
params.omegaA = 0.2;

% Choice of homophily function
% 'power'  -> F(h,L)=L+(1-L)h^alpha
% 'assort' -> F(h,L)=L+(1-L)h/(h+kappa(1-h))
params.homophily_type = 'power';
% params.homophily_type = 'assort';

%% Initial conditions
N0  = 10000;
SV0 = 5000;
SH0 = 3999;
SA0 = 1000;
V0  = 0;

IV0 = 0;
IH0 = 1;
IA0 = 0;
R0  = 0;

y0 = [SV0 SH0 SA0 V0 IV0 IH0 IA0 R0];

tspan = [0 500];

%% Sweep over homophily
hvals = linspace(0,1,100);
final_size = zeros(size(hvals));
final_vax  = zeros(size(hvals));

for k = 1:length(hvals)
    params.h = hvals(k);

    [t,y] = ode45(@(t,y) homophily_replicator_ode(t,y,params), tspan, y0);

    final_size(k) = 100*y(end,8)/N0;           % R(end)
    final_vax(k)  = 100*y(end,4)/N0;    % vaccinated percentage
end

%% Plot final vaccination coverage
 figure;
 hold on
 plot(hvals, final_vax, 'LineWidth', 3);
 xlabel('Homophily index h');
 ylabel('Final Vaccination Coverage (%)');
 grid on;
 
 set(gca, ...
     'Box', 'off', ...
     'TickDir', 'out', ...
     'LineWidth', 0.8, ...
     'FontSize', 14, ...
    'XColor', 'k', ...
     'YColor', 'k');
 set(gca,'xtick',0:0.1:1,'xticklabel',0:0.1:1);
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
% %% Plot final epidemic size
% figure;
hold on
plot(hvals, final_size, 'LineWidth', 3);
xlabel('Homophily index h');
ylabel('Attack Rate (%)');
grid on;

set(gca, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'LineWidth', 0.8, ...
    'FontSize', 14, ...
    'XColor', 'k', ...
    'YColor', 'k');
set(gca,'xtick',0:0.1:1,'xticklabel',0:0.1:1);
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%% Example: time dynamics for one value of h
% params.h = 0.1;
% [t,y] = ode45(@(t,y) homophily_replicator_ode(t,y,params), tspan, y0);
% 
% I_total = y(:,5) + y(:,6) + y(:,7);
% 
% % figure;
% hold on
% plot(t, I_total, 'LineWidth', 3);
% xlabel('Time');
% ylabel('Total infected');
% grid on;
% 
% set(gca, ...
%     'Box', 'off', ...
%     'TickDir', 'out', ...
%     'LineWidth', 0.8, ...
%     'FontSize', 14, ...
%     'XColor', 'k', ...
%     'YColor', 'k');