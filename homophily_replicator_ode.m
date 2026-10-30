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