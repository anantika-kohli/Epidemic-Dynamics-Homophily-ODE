clear; clc; close all;
% replicator + hill model in homophily_replicator_hill_ode.m.
% For both homophily forms (power, assort) and three seeds (IV, IH, IA):
% Plots vaccination coverage vs h
% Plots epidemic size (total infections) vs h

%% Base parameters 
P.beta  = 1.1e-5; P.c = 8; P.gamma = 1/5; P.xi = 0.025; P.eps = 0.7;
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

seedCols  = [5 6 7];                          % columns for IV, IH, IA
seedNames = {'I_V = 1','I_H = 1','I_A = 1'};
homTypes  = {'power','assort'};
homTitles = {'Power-law homophily','Assortativity-based homophily'};

%% Two scenarios: weak then strong adaptation 
% pv = rate of hesitant switching to pro-vax switching, pa = hesitant switchong to anti-vax.
pvSet = [0.05 0.30];     
paSet = [0.05 0.30];    
tagSet = {'weak adaptation','strong adaptation'};

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
                y0 = [5000 4000 1000 0 0 0 0 0];
                y0(seedCols(s))   = y0(seedCols(s))   + 1;                                    % seed one infected
                y0(seedCols(s)-4) = y0(seedCols(s)-4) - 1;                                    % remove from its susceptible pool
                [~,y] = ode45(@(t,y) homophily_replicator_hill_ode(t,y,P), tspan, y0);
                cov(k,s) = 100*y(end,4)/N0;                                                   % V(end): coverage (stock)
                siz(k,s) = y(end,8);                                                          % R(end): total infections
            end
        end

        figure(covFig); subplot(1,2,ht);
        plot(hvals, cov, 'LineWidth', 2.5);
        ylim([20 90]);
        xlabel('Homophily index h'); ylabel('Final Vaccination Coverage (%)');
        title(homTitles{ht}); legend(seedNames,'Location','best'); grid on; styleAxes;

        figure(sizFig); subplot(1,2,ht);
        plot(hvals, siz, 'LineWidth', 2.5);
        ylim([0 6000]);
        xlabel('Homophily index h'); ylabel('Total infections');
        title(homTitles{ht}); legend(seedNames,'Location','best'); grid on; styleAxes;
    end

    figure(covFig); sgtitle(['Vaccination coverage during an epidemic: ' tagSet{sc}]);
    figure(sizFig); sgtitle(['Effect of homophily on epidemic size: ' tagSet{sc}]);
end

%% Local function 
function styleAxes()
set(gca,'Box','off','TickDir','out','LineWidth',0.8,'FontSize',12, ...
    'XColor','k','YColor','k','XTick',0:0.2:1);
xlim([0 1]);
end
