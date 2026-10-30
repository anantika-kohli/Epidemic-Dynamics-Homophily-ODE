function [dydt, extras] = homophily_replicator_hill_ode(~, y, p)

% State order: [SV SH SA V IV IH IA R]

% State variables
SV = y(1); SH = y(2); SA = y(3); V = y(4);
IV = y(5); IH = y(6); IA = y(7); R = y(8);

% Homophily function F(h,L): fraction of within-group contacts
switch lower(p.homophily_type)
    case 'power'      % F = L + (1-L) h^alpha
        F = p.L + (1-p.L)*p.h^p.alpha;
    case 'assort'     % F = L + (1-L) h/(h + kappa(1-h))      
        F = p.L + (1-p.L)*( p.h / (p.h + p.kappa*(1-p.h) + eps) );
    case 'logistic'   % F = L + (1-L)/(1+exp(-k(h-h0)))
        F = p.L + (1-p.L)/(1+exp(-p.k*(p.h - p.h0)));
    otherwise
        error('homophily_type must be ''power'', ''assort'', or ''logistic''.');
end

% Forces of infection
lambdaV = p.beta*p.c*( F*IV + (1-F)*(IH+IA) );
lambdaH = p.beta*p.c*( F*IH + (1-F)*(IV+IA) );
lambdaA = p.beta*p.c*( F*IA + (1-F)*(IV+IH) );

% Prevalence-dependent risk perception omega_i(I)  
% omega_i(I) = omega_i_min + (omega_i_max - omega_i_min) * I^m / (I^m + I0^m)
I    = IV + IH + IA;                            % total infectious prevalence
hill = I^p.m / (I^p.m + p.I0^p.m + eps);        % Hill response in [0,1)
omegaV = p.omegaV_min + (p.omegaV_max - p.omegaV_min)*hill;
omegaH = p.omegaH_min + (p.omegaH_max - p.omegaH_min)*hill;
omegaA = p.omegaA_min + (p.omegaA_max - p.omegaA_min)*hill;

% Behavioural payoffs
PiV = p.bV - p.cV - omegaV*lambdaV;
PiH = p.bH - p.cH - omegaH*lambdaH;
PiA = p.bA - p.cA - omegaA*lambdaA;

% Attitude frequencies among the active population
M = SV + SH + SA + V + IV + IH + IA;
if M > 0
    xV = (SV + V + IV)/M;
    xA = (SA + IA)/M;
else
    xV = 0; xA = 0;
end

% Replicator switching 
THV = p.pv*p.c*(1-F)*SH*xV*max(PiV - PiH, 0);
THA = p.pa*p.c*(1-F)*SH*xA*max(PiA - PiH, 0);

% Full dynamical system 
dSV = -lambdaV*SV + THV - p.xi*SV;
dSH = -lambdaH*SH - THV - THA;
dSA = -lambdaA*SA + THA;
dV  =  p.xi*SV - (1-p.eps)*lambdaV*V;
dIV =  lambdaV*SV + (1-p.eps)*lambdaV*V - p.gamma*IV;
dIH =  lambdaH*SH - p.gamma*IH;
dIA =  lambdaA*SA - p.gamma*IA;
dR  =  p.gamma*(IV+IH+IA);

dydt = [dSV; dSH; dSA; dV; dIV; dIH; dIA; dR];
if nargout>1
    extras.PiV = PiV;  extras.PiH = PiH;  extras.PiA = PiA;
    extras.THV = THV;  extras.THA = THA;
    extras.F   = F;
    extras.lambdaV = lambdaV; extras.lambdaH = lambdaH; extras.lambdaA = lambdaA;
end
end

