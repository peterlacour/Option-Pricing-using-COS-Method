function phi_fmls = fmls_char_fn(alpha, omega, eta, T, mu)

%{
 This code computes the Characteristic Function for the Finite Moment Log
 Stable Levy Process. Notation follows Carr and Wu (2003), eq. 32, p. 8

 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (21.03.2019)

 phi_fmls = fmls_char_fn(alpha, omega, sigma, T)


 Inputs : alpha         - tail index, element of (1,2]
        : eta           - dispersion parameter
        : mu            - r - q
        : T             - time to maturity
        : omega         - argument of the characteristic function


Outputs : phi_fmls      - characteristic function value

%}

w = eta^alpha * sec(pi*alpha/2);

% Page 29 equation (5.11)
phi_fmls    = exp(1i * omega * (mu + w) * T - ( 1i * omega * eta ).^alpha * T * sec(pi*alpha/2));

