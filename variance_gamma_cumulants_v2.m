function [c1, c2, c4, w] = variance_gamma_cumulants_v2( u_0, T, theta, eta, mu )

%{
 This code computes up to the 2nd cumulant of ln(St/K) for the Variance Gamma Model
 Equations are given in Fang and Oosterlee (2008), Table 11, p. 21

 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (21.03.2019)

 [c1, c2, omega] = vg_cumulants_v1( u_0, T, theta, eta, mu )


 Inputs : u_0           - variance
        : T             - time to maturity
        : theta         - drift term
        : mu            - mean return
        : eta           - variance rate

Outputs : c1            - first cumulant (mean)
        : c2            - second cumulant (variance)
        : omega         - drift correction term

%}


% First cumulant (mean)
c1 = ( mu + theta ) * T;

% Second cumulant (variance)
c2 = ( u_0 + eta * theta^2 ) * T;

% Second cumulant (variance)
c4 = 3 * ( u_0^2 * eta + 2 * theta^4 * eta^3 + 4 * u_0 * theta^2 * eta^2 ) * T;

% Drift correction term (0 for the Heston Model)
w = (1 / eta) * log( (1 - theta * eta - u_0*eta) / 2 );