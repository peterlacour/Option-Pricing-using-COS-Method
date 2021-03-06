function phi_cgmy = cgmy_char_fn(mu, u_0, C, G, Y, M, a, b, k, T)

%{

 This code computes the Characteristic Function for the CGMY model process. 
 Notation follows Fang and Oosterlee (2008) eq. 31, p. 8 and parameters are
 inspired by Carr et al. (2002).

 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (22.03.2019)
         : 2.0 (23.03.2019) Added internal computation of omegas

 phi_cgmy = cgmy_char_fn(mu, u_0, C, G, Y, M, a, b, k, T)


 Inputs : mu            - r - q
        : T             - time to maturity
        : omega         - argument of the characteristic function
        : C             - measure of the overall level of activity of the
                          process, controls the kurtosis
        : G             - controls the exponential decay on the right tail 
                          of the process
        : M             - controls the exponential decay on the left tail 
                          of the process
        : Y             - characterizes the fine structure of the jumps of 
                          the CGMY process
        : u_0           - variance

Outputs : phi_cgmy      - characteristic function values [0:N-1] vector

%}


% Vector of N evaluation arguments for the characteristic function
omega   = k .* pi / (b - a);


% Characteristic function
phi_cgmy = exp( 1i * omega * ( mu ) * T - 0.5 * omega.^2 * u_0 * T) .* ...
           exp( T * C * gamma( -Y ) * ( (M - 1i * omega).^Y - M^Y + (G + 1i * omega).^Y - G^Y ) );



