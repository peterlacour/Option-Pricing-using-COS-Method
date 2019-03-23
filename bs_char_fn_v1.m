function phi_bs = bs_char_fn_v1(mu, u_0, a, b, k, T)

%{
 This code computes the Characteristic Function for the Heston Model
 Notation follows Fang and Oosterlee (2008), eq. 32, p. 8

 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (21.03.2019)
         : 2.0 (23.03.2019) Added internal computation of omegas

 phi_hest = heston_char_fn_v2(mu, lambda, u_bar, u_0, eta, rho, a, b, k, T)


 Inputs : mu            - log price drift rate
        : lambda        - speed of mean reversion
        : u_bar         - mean (long run) volatility
        : u_0           - initial volatility
        : eta           - volatility of the volatility (vol of vol)
        : rho           - correlation between Wiener processes (W1 and W2)
        : T             - time to maturity
        : a             - Cosine argument (lower truncation bound) 
        : b             - Cosine argument (upper truncation bound)
        : k             - Vector of N evaluation intervals

Outputs : phi_hest      - characteristic function values [0:N-1] vector

%}


% Vector of N evaluation arguments for the characteristic function
omega   = k .* pi / (b - a);

phi_bs = exp((mu - 0.5 * u_0) * 1i * T .* omega - 0.5 * u_0 * T * omega .^2);
