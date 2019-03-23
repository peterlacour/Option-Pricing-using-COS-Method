%% Housekeeping

clear, clc;
%% Model Parameters

% Heston Model Parameters [defaults from Fang and Oosterlee (2008)]
lambda      =  1.5768;                          % Mean Reversion Speed
u_bar       =  0.0398;                          % Long-run Variance
u_0         =  0.0175;                          % Initial Variance
eta         =  0.5751;                          % Volatility of the Volatility Parameter
rho         = -0.5711;                          % Correlation between Wiener Processes


% VG Parameters
theta = -0.1436;


% CGMY Parameters from CGMY (2003)
C = 24.79;
G = 94.45;
Y = 0.2495;
M = 95.79;


% Option / Underlying / Market Parameters
T    = 1;                                       % Time to Maturity
r    = 0;                                       % Risk-free Rate
q    = 0;                                       % Dividend Yield
S0   = 100;                                     % Initial Underlying Price
mu   = r - q;                                   % Price Drift Rate


%% COS-FFT Truncation Bounds

% Heston cumulants and integration bounds
[c1, c2, ~]     = heston_cumulants_v1(mu, lambda, u_bar, u_0, eta, rho, T);
[a_hest, b_hest]= cos_truncation_range_v2(c1,c2,0,12);

% CGMY cumulants and integration bounds
[c1, c2, c4, ~] = cgmy_cumulants_v2( u_0, T, mu, C, G, M, Y);
[a_cgmy, b_cgmy]= cos_truncation_range_v2(c1,c2,c4,10);

% VG cumulants and integration bounds
[c1, c2, c4, ~] = variance_gamma_cumulants_v2( u_0, T, theta, 1, mu );
[a_vg, b_vg]    = cos_truncation_range_v2(c1,c2,c4,10);

% BS cumulants and integration bounds
[c1, c2, c4, ~] = bs_cumulants_v1(u_0, mu, T );
[a_bs, b_bs]    = cos_truncation_range_v2(c1,c2,c4,10);


%% COS - FFT Parameters

N       = 5000;                 % Number of points to evaluate
k       = 0:(N - 1);            % Vector of N evaluation intervals
K       = 70:130;               % Vector of M strike prices to evaluate
x       = log(S0 ./ K);         % Vector of M log prices


% Characteristic functions for the log stock price
phi_hest     = heston_char_fn_v2(mu, lambda, u_bar, u_0, eta, rho, a_hest, b_hest, k, T);
phi_cgmy     = cgmy_char_fn(mu, u_0, C, G, Y, M, a_cgmy, b_cgmy, k, T);
phi_vg       = vg_char_fn(u_0, theta, a_vg, b_vg, k, T);
phi_bs       = bs_char_fn_v1(mu, u_0, a_bs, b_bs, k, T);


% FFT-COS Prices
[C_COS_hest, P_COS_hest] = cos_option_price_v1(a_hest, b_hest, k, K, phi_hest, x, mu, T);
[C_COS_cgmy, P_COS_cgmy] = cos_option_price_v1(a_cgmy, b_cgmy, k, K, phi_cgmy, x, mu, T);
[C_COS_vg,   P_COS_vg]   = cos_option_price_v1(a_vg,   b_vg,   k, K, phi_vg,   x, mu, T);
[C_COS_bs,   P_COS_bs]   = cos_option_price_v1(a_bs,   b_bs,   k, K, phi_bs,   x, mu, T);



%% Plots

% Black and Scholes Option Prices
% Based on code by Peter.Gruber@unisg.ch, and Paul.Soderlind@unisg.ch
[C_BS, P_BS, ~, ~] = black_scholes_price(S0,K,r,T,sqrt(u_0),q);

subplot(2,1,1)
plot(K,C_COS_hest,'r'), grid on, hold on;
plot(K,C_BS,'--k'), hold off;

subplot(2,1,2)
plot(K,C_COS_hest-C_BS'), grid on;


subplot(2,1,1)
plot(K,P_COS_hest,'r'), grid on, hold on;
plot(K,P_BS,'--k'), hold off;

subplot(2,1,2)
plot(K,P_COS_hest-P_BS'), grid on;


%% SETUP
Npoints = 5000;
N = 5000;                         % default value 25, try 100,500,5000

% Heston Plot
% Heston truncation bounds
figure
MakePdfPlot2(phi_hest,a_hest,b_hest,k,mu,u_0,N,1, 'Heston Pdf');


% Variance Gamma Plot
MakePdfPlot2(phi_vg,a_vg,b_vg,k,mu,u_0,N,2, 'Variance-Gamma Pdf');


% CGMY Plot
% CGMY - Parameters from CGMY ( 2003 )
MakePdfPlot2(phi_cgmy,a_cgmy,b_cgmy,k,mu,u_0,N,3, 'CGMY Pdf');





%% ========== PRIVATE FUNCTIONS ==========
function MakePdfPlot2(cf,a,b,k,mu,u_0,N,plotpos,name)
bma = b-a;
x = linspace(a,b,N);
true=normpdf(x,mu, sqrt(u_0) );

Fk = 2/bma * real( cf.*exp(-1i*k*a*pi/bma)   );
Fk(1)=0.5*Fk(1);

for l=1:length(x)
    pdf(l)=sum(Fk.*cos (k*pi*(x(l)-a)/bma) );
end

subplot(2,3,plotpos)
hold on
plot(x,pdf, 'DisplayName', name), grid on;
title(name)
plot(x, true, 'DisplayName', 'Normal Pdf')
legend
hold off

subplot(2,3,plotpos+3)
plot(x,pdf(:)-true(:)), grid on;
title(strcat(name, " - Normal Pdf"))

end
