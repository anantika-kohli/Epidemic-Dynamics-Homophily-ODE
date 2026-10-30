clear; clc; close all;
% ==================================================================
%  Reproduces RUN8-RUN12 from your sensitivity log.
%
%  Model = pv/pa direction-specific switching rates (as in
%  homophily_replicator_hill_ode.m) BUT with CONSTANT omega_i values
%  (as in homophily_replicator_ode.m / trial1.m), not the Hill function.
%  This matches the model version that was actually in use when
%  RUN8-RUN12 were generated (pv/pa already existed by RUN4-6, but the
%  Hill-type omega_i(I) was a later addition on top of the "original"
%  constant-omega code).
%
%  RUN_ID below applies the one-parameter change for that run.
%  Set RUN_ID = 0 to run the plain baseline (no override).
% ==================================================================

RUN_ID = 12;   % 0 = baseline, or 8 / 9 / 10 / 11 / 12

%% ---------------- Base parameters (confirmed from your files) ----------------
P.beta  = 1.1e-5; P.c = 8; P.gamma = 1/5; P.xi = 0.025; P.eps = 0.7;
P.L = 0.7; P.alpha = 2; P.kappa = 0.2;
P.bV = 0.20; P.bH = 0.10; P.bA = 0.15;
P.cV = 0.08; P.cH = 0.03; P.cA = 0.02;

% Constant omega values (from trial1.m -- the pre-Hill original code)
P.omegaV = 1.5;
P.omegaH = 0.8;
P.omegaA = 0.2;

N0    = 10000;
tspan = [0 500];
hvals = linspace(0,1,80);

% --- Apply the one-parameter change for the selected RUN ---
switch RUN_ID
    case 8,  P.bA = 0.20;   % RUN8:  increase benefit for anti (from 0.15)
    case 9,  P.bH = 0.16;   % RUN9:  increase benefit for hesitant (from 0.10)
    case 10, P.cH = 0.06;   % RUN10: increase cost for hesitant (from 0.03)
    case 11, P.cA = 0.05;   % RUN11: increase cost for anti (from 0.02)
    case 12, P.bH = 0.05;   % RUN12: decrease benefit for hesitant (from 0.10)
    otherwise                % RUN_ID = 0 -> plain baseline, no override
end

% --- pv/pa weak/strong (confirmed default from MainRun_replicator_figures_1.m) ---
pvSet = [0.05 0.30];   % weak, strong -> pro-vax switching rate
paSet = [0.05 0.30];   % weak, strong -> anti-vax switching rate
tagSet = {'weak adaptation','strong adaptation'};

% --- Seeding: exact base split from your MainRun/trial1 files ---
seedCols  = [5 6 7];                 % IV, IH, IA columns in y
seedNames = {'I_V = 1','I_H = 1','I_A = 1'};
homTypes  = {'power','assort'};
homTitles = {'Power-law homophily','Assortativity-based homophily'};

%% ---------------- Two scenarios: weak then strong adaptation ----------------
for sc = 1:2
    P.pv = pvSet(sc);
    P.pa = paSet(sc);
    covFig = figure('Name',['Coverage - ' tagSet{sc}],'Position',[100 100 1100 430]);
    sizFig = figure('Name',['Epidemic size - ' tagSet{sc}],'Position',[100 100 1100 430]);

    for ht = 1:2
        P.homophily_type = homTypes{ht};
        cov = zeros(numel(hvals),3);
        siz = zeros(numel(hvals),3);

        for s = 1:3
            for k = 1:numel(hvals)
                P.h = hvals(k);
                y0 = [5000 4000 1000 0 0 0 0 0];      % SV SH SA V IV IH IA R
                y0(seedCols(s))   = y0(seedCols(s))   + 1;  % seed one infected
                y0(seedCols(s)-4) = y0(seedCols(s)-4) - 1;  % remove from its susceptible pool

                [~,y] = ode45(@(t,y) homophily_replicator_pvpa_const(t,y,P), tspan, y0);

                cov(k,s) = 100*y(end,4)/N0;   % V(end): coverage (stock)
                siz(k,s) = y(end,8);          % R(end): total infections
            end
        end

        figure(covFig); subplot(1,2,ht);
        plot(hvals, cov, 'LineWidth', 2.5);
        xlabel('Homophily index h'); ylabel('Final Vaccination Coverage (%)');
        title(homTitles{ht}); legend(seedNames,'Location','best'); grid on; styleAxes();

        figure(sizFig); subplot(1,2,ht);
        plot(hvals, siz, 'LineWidth', 2.5);
        xlabel('Homophily index h'); ylabel('Total infections');
        title(homTitles{ht}); legend(seedNames,'Location','best'); grid on; styleAxes();
    end

    figure(covFig); sgtitle(['Vaccination coverage during an epidemic: ' tagSet{sc}]);
    figure(sizFig); sgtitle(['Effect of homophily on epidemic size: ' tagSet{sc}]);
end

%% ---------------- Local functions ----------------
function styleAxes()
set(gca,'Box','off','TickDir','out','LineWidth',0.8,'FontSize',12, ...
    'XColor','k','YColor','k','XTick',0:0.2:1);
xlim([0 1]);
end

function dydt = homophily_replicator_pvpa_const(~, y, p)
% Same structure as homophily_replicator_ode.m (constant omega_i), but
% with direction-specific pv/pa switching rates instead of one shared rho.
SV = y(1); SH = y(2); SA = y(3); V = y(4);
IV = y(5); IH = y(6); IA = y(7); %#ok<NASGU>

switch lower(p.homophily_type)
    case 'power'
        F = p.L + (1-p.L)*p.h^p.alpha;
    case 'assort'
        F = p.L + (1-p.L)*( p.h / (p.h + p.kappa*(1-p.h) + eps) );
    otherwise
        error('Unknown homophily_type. Use ''power'' or ''assort''.');
end

lambdaV = p.beta*p.c*( F*IV + (1-F)*(IH+IA) );
lambdaH = p.beta*p.c*( F*IH + (1-F)*(IV+IA) );
lambdaA = p.beta*p.c*( F*IA + (1-F)*(IV+IH) );

M = SV+SH+SA+V+IV+IH+IA;
if M > 0
    xV = (SV+V+IV)/M;
    xA = (SA+IA)/M;
else
    xV = 0; xA = 0;
end

PiV = p.bV - p.cV - p.omegaV*lambdaV;
PiH = p.bH - p.cH - p.omegaH*lambdaH;
PiA = p.bA - p.cA - p.omegaA*lambdaA;

THV = p.pv*p.c*(1-F)*SH*xV*max(PiV-PiH, 0);
THA = p.pa*p.c*(1-F)*SH*xA*max(PiA-PiH, 0);

dSV = -lambdaV*SV + THV - p.xi*SV;
dSH = -lambdaH*SH - THV - THA;
dSA = -lambdaA*SA + THA;
dV  =  p.xi*SV - (1-p.eps)*lambdaV*V;
dIV =  lambdaV*SV + (1-p.eps)*lambdaV*V - p.gamma*IV;
dIH =  lambdaH*SH - p.gamma*IH;
dIA =  lambdaA*SA - p.gamma*IA;
dR  =  p.gamma*(IV+IH+IA);

dydt = [dSV; dSH; dSA; dV; dIV; dIH; dIA; dR];
end
