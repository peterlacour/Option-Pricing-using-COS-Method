function phi_cgmy = cgmy_char_fn(mu, u_0, C, G, Y, M, a, b, k, T)

%{

 This code computes the Characteristic Function for the Finite Moment Log
 Stable Levy Process. Notation follows Carr and Wu (2003), eq. 32, p. 8

 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (22.03.2019)
         : 2.0 (23.03.2019) Added internal computation of omegas

 phi_cgmy = cgmy_char_fn(mu, u_0, C, G, Y, M, a, b, k, T)


 Inputs : mu            - r - q
        : T             - time to maturity
        : omega         - argument of the characteristic function
        : C             - Carr Parameter
        : G             - Geman Parameter
        : M             - Madan Parameter
        : Y             - Yor Parameter
        : u_0           - variance

Outputs : phi_cgmy      - characteristic function values [0:N-1] vector

%}


% Vector of N evaluation arguments for the characteristic function
omega   = k .* pi / (b - a);


% Characteristic function
phi_cgmy = exp( 1i * omega * ( mu ) * T - 0.5 * omega.^2 * u_0 * T) .* ...
           exp( T * C * gamma( -Y ) * ( (M - 1i * omega).^Y - M^Y + (G + 1i * omega).^Y - G^Y ) );



