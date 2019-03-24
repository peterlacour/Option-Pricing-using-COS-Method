function phi_bs = bs_char_fn_v1(u_0, a, b, k, T)

%{
 This code computes the Characteristic Function for the Black-Scholes Model
 Notation follows Schmelzle (2010), eq. 5.10, p. 37

 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (21.03.2019)
         : 2.0 (23.03.2019) Added internal computation of omegas

 phi_bs = bs_char_fn_v1(u_0, a, b, k, T)

 Inputs : u_0           - variance
        : T             - time to maturity
        : a             - Cosine argument (lower truncation bound) 
        : b             - Cosine argument (upper truncation bound)
        : k             - Vector of N evaluation intervals

Outputs : phi_bs        - characteristic function values [0:N-1] vector

%}


% Vector of N evaluation arguments for the characteristic function
omega   = k .* pi / (b - a);

phi_bs = exp(( -0.5 * u_0 ) * 1i * T .* omega - 0.5 * u_0 * T * omega .^2);
