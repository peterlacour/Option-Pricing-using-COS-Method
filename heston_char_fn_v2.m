function phi_hest = heston_char_fn_v2(mu, lambda, u_bar, u_0, eta, rho, a, b, k, T)

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


% D and G parameters: Fang and Oosterlee (2008) p. 8
D = sqrt( (lambda - 1i * rho * eta * omega) .^ 2 + (omega .^ 2 + 1i * omega) * eta ^ 2 );
G = (lambda - 1i * rho * eta * omega - D) ./ (lambda - 1i * rho * eta * omega + D);


% Characteristic function for the Heston Model: 
% Fang and Oosterlee (2008) p. 8
phi_hest = exp( 1i * omega * mu * T + ...
           (u_0 / eta ^ 2) * (1 - exp(-D * T)) ./ (1 - G .* exp(-D * T)) .* (lambda - 1i * rho * eta * omega - D) + ...
           (lambda * u_bar / eta ^ 2 * (T * (lambda - 1i * rho * eta * omega - D) - 2 * log( (1 - G .* exp(-D * T)) ./ (1 - G) ))));
      
  