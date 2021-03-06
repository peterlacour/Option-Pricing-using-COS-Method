function [C_COS, P_COS] = cos_option_price_v1(a, b, k, K, phi, x, mu, T)
%{
 This code computes the call or put option price using the COS-FFT Method
 described in Fang and Oosterlee (2008).

 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (21.03.2019)

 [c1, c2, omega] = heston_cumulants(mu, lambda, u_bar, u_0, eta, rho, T)


 Inputs : a             - lower truncation bound
        : b             - upper truncation bound
        : k             - Vector of N evaluation intervals
        : K             - Vector of M strike prices to evaluate
        : phi           - values from given characteristic function
        : x             - Vector of M log prices 
        : mu            - drift = r - q 
        : T             - time to maturity


Outputs : C_COS         - call option price
        : P_COS         - put option price


% 
K       = strikes';             % Vector of M strike prices to evaluate
%}

% Vk Integrals vector [0:N-1]
[V_k_call, V_k_put] = cos_series_analytical_integrals_v1(a, b, k, K);


% (Vectorised) Approximate cosine expansion with approximated characteristic function
% Fang and Oosterlee (2008) eq. 9, page 4
F_k     = real((ones(size(K, 2), 1) * phi) .* exp(1i .* (x - a)' * k .* pi ./ (b - a))); 


% Weigh First term 1/2
F_k(:,1)  = 0.5 * F_k(:,1);


% FFT-COS Prices
% Fang and Oosterlee (2008), eq. 24, page 6
C_COS = sum(F_k .* V_k_call, 2) * exp(-mu * T);
% Fang and Oosterlee (2008), eq. 25, page 6
P_COS = sum(F_k .* V_k_put, 2)  * exp(-mu * T);