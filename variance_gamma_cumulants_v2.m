function [c1, c2, c4, w] = variance_gamma_cumulants_v2( u_0, T, theta, mu )

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


Outputs : c1            - first cumulant (mean)
        : c2            - second cumulant (variance)
        : omega         - drift correction term

%}
% variance rate set to 1 for now
v = 1

% First cumulant (mean)
c1 = ( mu + theta ) * T;

% Second cumulant (variance)
c2 = ( u_0 + v * theta^2 ) * T;

% Second cumulant (variance)
c4 = 3 * ( u_0^2 * v + 2 * theta^4 * v^3 + 4 * u_0 * theta^2 * v^2 ) * T;

% Drift correction term (0 for the Heston Model)
w = (1 / v) * log( (1 - theta * v - u_0*v) / 2 );