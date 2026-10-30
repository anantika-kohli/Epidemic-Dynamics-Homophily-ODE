function [cov, siz] = run_scenario(P, homophily_type, pv, pa, hvals, seedCols, tspan, N0)
% RUN_SCENARIO  Core simulation engine, reused by every figure.
%
% Sweeps h across hvals for one (homophily_type, pv, pa) combination and
% returns:
%   cov(k,s) = final vaccination coverage (%) for seed s at hvals(k)
%   siz(k,s) = total infections (cumulative R at t_end) for seed s at hvals(k)
%
% Uses homophily_replicator_hill_ode.m unchanged -- keep all three files
% in the same folder.

P.homophily_type = homophily_type;
P.pv = pv;
P.pa = pa;

nS = numel(seedCols);
cov = zeros(numel(hvals), nS);
siz = zeros(numel(hvals), nS);

for s = 1:nS
    for k = 1:numel(hvals)
        P.h = hvals(k);
        y0 = [5000 4000 1000 0 0 0 0 0];
        y0(seedCols(s))   = y0(seedCols(s))   + 1;   % seed one infected
        y0(seedCols(s)-4) = y0(seedCols(s)-4) - 1;   % remove from its susceptible pool
        [~, y] = ode45(@(t,y) homophily_replicator_hill_ode(t,y,P), tspan, y0);
        cov(k,s) = 100*y(end,4)/N0;   % V(end): final coverage
        siz(k,s) = y(end,8);          % R(end): total infections
    end
end
end
