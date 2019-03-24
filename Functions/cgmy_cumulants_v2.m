function [c1, c2, c4, w] = cgmy_cumulants_v2( u_0, T, mu, C, G, M, Y)

%{
 This code computes the first, second and fourth cumulant of ln(St/K) 
 for the Variance Gamma Model Equations are given in Fang and Oosterlee (2008), Table 11, p. 21

 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (21.03.2019)

 [c1, c2, omega] = vg_cumulants_v1( u_0, T, theta, eta, mu )


 Inputs : u_0           - variance
        : T             - time to maturity
        : mu            - mean return
        : eta           - variance rate

Outputs : c1            - first cumulant (mean)
        : c2            - second cumulant (variance)
        : w             - drift correction term

%}


% First cumulant (mean)
c1 = mu * T + C * T * gamma(1-Y) * (M^(Y-1) - G^(Y-1));

% Second cumulant (variance)
c2 = u_0 * T + C * T * gamma(2-Y) * (M^(Y-2) - G^(Y-2));

% Fourth cumulant
c4 = C * T * gamma(4-Y) * ( M^(Y-4) + G^(Y-4) );
 
 % Drift correction term (0 for the Heston Model)
 w = - C * gamma(-Y) * (( M-1)^Y - M^Y + (G+1)^Y - G^Y );
 