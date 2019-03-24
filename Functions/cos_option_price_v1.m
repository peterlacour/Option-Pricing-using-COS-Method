function [C_COS, P_COS] = cos_option_price_v1(a, b, k, K, phi, x, mu, T)


% Vk Integrals vector [0:N-1]
[V_k_call, V_k_put] = cos_series_analytical_integrals_v1(a, b, k, K);


% (Vectorised) Approximate cosine expansion with approximated characteristic function
F_k     = real((ones(size(K, 2), 1) * phi) .* exp(1i .* (x - a)' * k .* pi ./ (b - a)));


% Weigh First term 1/2
F_k(:,1)  = 0.5 * F_k(:,1);


% FFT-COS Prices
C_COS = sum(F_k .* V_k_call, 2) * exp(-mu * T);
P_COS = sum(F_k .* V_k_put, 2)  * exp(-mu * T);