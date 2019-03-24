function phi_vg = vg_char_fn(u_0, theta, a, b, k, T, v, S0, mu)

%{
 This code computes the Characteristic Function for the Variance Gamma Process. 
 Notation follows Madan, Carr and Chang (1998), eq. 7, p. 8. 
 theta's average value in Madan, Carr, Chang (1998) is -0.1436, 
 theta = 0 gives the symmetric variance gamma,
 theta < 0 incorporates negatively skewness in the model.
 
 Source: A.Itkin ”Pricing options with VG model using FFT”. 
                   The Variance Gamma and Related Financial Models. August 9, 2007 – p. 7

 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (21.03.2019)
         : 2.0 (23.03.2019) Added internal computation of omegas

 phi_vg = vg_char_fn(u_0, v, theta, a, b, k, T)


 Inputs : u_0           - variance
        : theta         - drift term
        : T             - time to maturity
        : omega         - argument of the characteristic function
        : v             - variance rate
        : a             - Cosine argument (lower truncation bound) 
        : b             - Cosine argument (upper truncation bound)
        : k             - Vector of N evaluation intervals


Outputs : phi_vg        - characteristic function values [0:N-1] vector

%}

% Vector of N evaluation arguments for the characteristic function
omega   = k .* pi / (b - a);

w = ( 1 / v ) * log(1 - theta * v - 0.5 * u_0 * v );

phi_vg = exp( log(S0 + ( mu + w ) * T ) ./ ( 1 - 1i * theta * v * omega + ( u_0 * v / 2 ) * omega.^2 ).^(T/v) );



