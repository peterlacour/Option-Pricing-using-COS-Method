function phi_vg = vg_char_fn(u_0, eta, theta, a, b, k, T)

%{
 This code computes the Characteristic Function for the Variance Gamma Process. 
 Notation follows Madan, Carr and Chang (1998), eq. 7, p. 8. 
 theta's average value in Madan, Carr, Chang (1998) is -0.1436, 
 theta = 0 gives the symmetric variance gamma,
 theta < 0 incorporates negatively skewness in the model.


 Authors : Baldi Lanfranchi, Federico
         : La Cour, Peter

 Version : 1.0 (21.03.2019)
         : 2.0 (23.03.2019) Added internal computation of omegas

 phi_vg = vg_char_fn(u_0, eta, theta, a, b, k, T)


 Inputs : u_0         - variance
        : theta         - drift term
        : T             - time to maturity
        : omega         - argument of the characteristic function
        : eta           - variance rate
        : a             - Cosine argument (lower truncation bound) 
        : b             - Cosine argument (upper truncation bound)
        : k             - Vector of N evaluation intervals


Outputs : phi_vg        - characteristic function values [0:N-1] vector

%}

% Vector of N evaluation arguments for the characteristic function
omega   = k .* pi / (b - a);

phi_vg = ( 1 ./ ( 1 - 1i * theta * eta * omega + ( u_0 * eta / 2 ) * omega.^2 ) ).^(T/eta);



