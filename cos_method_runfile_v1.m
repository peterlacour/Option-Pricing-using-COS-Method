%{
 
 Title         : COS-FFT Method to Option Pricing using different
                 Characteristic Functions

 Authors       : Baldi Lanfranchi, Federico
               : La Cour, Peter
 Version       : 1.0 (24.03.2019)

 Place, Time   : St. Gallen, 24.03.2019
 Description   : The code includes five characteristic functions: the
                 Black-Scholes, the Heston, the Variance Gamma, the CGMY 
                 and the FMLS (unused) and compares the distributional 
                 characteristics by computing the distributions to a normal 
                 distribution using the cosine series Fourier expansion . 
                 The code then prices an option chain of the SPX with 
                 maturity date of 20.03.2020 using the different 
                 characteristic functions and compares it to the analytical 
                 solution of the Black-Scholes model.

Future       
Improvements  :  Optimise parameters and compare pricing method to real
                 option prices
%}

%% Housekeeping
clear, clc;
addpath("./Functions")

%% Import Option Chain Data 

parameters = readtable("params.csv");
strikes    = readtable("strike_prices.csv");
optionBids = readtable("option_prices.csv");

parameters = table2array(parameters(1,:));
strikes    = table2array(strikes(:,1));
optionBids = table2array(optionBids(:,1));

%% Model Parameters

% Heston Model Parameters [defaults from Fang and Oosterlee (2008)]
lambda      =  1.5768;                          % Mean Reversion Speed
u_bar       =  0.0398;                          % Long-run Variance
u_0         =  0.0175;                          % Initial Variance
eta         =  0.5751;                          % Volatility of the Volatility Parameter
rho         = -0.5711;                          % Correlation between Wiener Processes


% VG Parameters [defaults from Madan et al. (1998) and Carr et al. (2002)]
theta      = -0.1436;
v          = 0.0403;


% CGMY Parameters [defaults from Carr et al. (2002)]
C          = 24.79;
G          = 94.45;
Y          = 0.2495;
M          = 95.79;

% Option / Underlying / Market Parameters
S0         = parameters(1);                           % Initial Underlying Price
r          = parameters(2);                           % Risk-free Rate
q          = parameters(3);                           % Dividend Yield
T          = parameters(4);                           % Time to Maturity
mu         = r - q;                                   % Price Drift Rate



%% COS-FFT Truncation Bounds
% Start timer
tic

% Heston cumulants and integration bounds
[c1, c2, ~]     = heston_cumulants_v1(mu, lambda, u_bar, u_0, eta, rho, T);
[a_hest, b_hest]= cos_truncation_range_v2(c1,c2,0,12);

% CGMY cumulants and integration bounds
[c1, c2, c4, ~] = cgmy_cumulants_v2( u_0, T, mu, C, G, M, Y);
[a_cgmy, b_cgmy]= cos_truncation_range_v2(c1,c2,c4,10);

% VG cumulants and integration bounds
[c1, c2, c4, ~] = variance_gamma_cumulants_v2( u_0, T, theta, mu, v );
[a_vg, b_vg]    = cos_truncation_range_v2(c1,c2,c4,10);

% BS cumulants and integration bounds
[c1, c2, c4, ~] = bs_cumulants_v1(u_0, mu, T );
[a_bs, b_bs]    = cos_truncation_range_v2(c1,c2,c4,10);


%% COS - FFT Parameters

N             = 500;                 % Number of points to evaluate
k             = 0:(N - 1);            % Vector of N evaluation intervals
K             = strikes';             % Vector of M strike prices to evaluate
x             = log(S0 ./ K);         % Vector of M log prices


% Characteristic functions for the log stock price
phi_hest     = heston_char_fn_v2(mu, lambda, u_bar, u_0, eta, rho, a_hest, b_hest, k, T);
phi_cgmy     = cgmy_char_fn(mu, u_0, C, G, Y, M, a_cgmy, b_cgmy, k, T);
phi_vg       = vg_char_fn(u_0, theta, a_vg, b_vg, k, T, v, mu);
phi_bs       = bs_char_fn_v1(u_0, a_bs, b_bs, k, T);


% FFT-COS Prices
[C_COS_hest, P_COS_hest] = cos_option_price_v1(a_hest, b_hest, k, K, phi_hest, x, mu, T);
[C_COS_cgmy, P_COS_cgmy] = cos_option_price_v1(a_cgmy, b_cgmy, k, K, phi_cgmy, x, mu, T);
[C_COS_vg,   P_COS_vg]   = cos_option_price_v1(a_vg,   b_vg,   k, K, phi_vg,   x, mu, T);
[C_COS_bs,   P_COS_bs]   = cos_option_price_v1(a_bs,   b_bs,   k, K, phi_bs,   x, mu, T);
% Stop timer
toc

% Black and Scholes Analytical Option Prices
% Based on code by Peter.Gruber@unisg.ch, and Paul.Soderlind@unisg.ch
[C_BS, P_BS, ~, ~] = black_scholes_price(S0,K,r,T,sqrt(u_0),q);


%% Distribution Plots
Npoints = 5000;
N       = 5000;                         

% BS Plot
figure
MakePdfPlot2(phi_bs,a_bs,b_bs,k,mu,u_0,T,N,1, 'B-S Pdf');

% Heston Plot
% Heston truncation bounds
MakePdfPlot2(phi_hest,a_hest,b_hest,k,mu,u_0,T,N,2, 'Heston Pdf');

% Variance Gamma Plot
MakePdfPlot2(phi_vg,a_vg,b_vg,k,mu,u_0,T,N,3, 'Variance-Gamma Pdf');

% CGMY Plot
% CGMY - Parameters from CGMY ( 2003 )
MakePdfPlot2(phi_cgmy,a_cgmy,b_cgmy,k,mu,u_0,T,N,4, 'CGMY Pdf');

%% Option Price Plots

% Call Plots
% BS
subplot(2,4,1)
plot(K,C_COS_bs,'r', 'DisplayName', 'BS'), grid on, hold on;
plot(K,C_BS,'--k', 'DisplayName', 'Analytical BS'),
title('Black-Scholes Model')
legend
axis([min(strikes) max(strikes) 0.0 1800])
subplot(2,4,5)
plot(K,(C_COS_bs-C_BS')), grid on, hold on;
title("Diff. between Cosine Method Black-Scholes" + newline + "and Analytical Black-Scholes")
axis([min(strikes) max(strikes) -5.0 70.0])

% Heston
subplot(2,4,2)
plot(K,C_COS_hest,'r', 'DisplayName', 'Heston'), grid on, hold on;
plot(K,C_BS,'--k', 'DisplayName', 'Analytical BS'),
title('Heston Model')
axis([min(strikes) max(strikes) 0.0 1800])
legend
subplot(2,4,6)
plot(K,(C_COS_hest-C_BS')), grid on, hold on;
title("Diff. between Heston" + newline + "and Analytical Black-Scholes")
axis([min(strikes) max(strikes) -5.0 70.0])

% Variance Gamma
subplot(2,4,3)
plot(K,C_COS_vg,'r', 'DisplayName', 'VG'), grid on, hold on;
plot(K,C_BS,'--k', 'DisplayName', 'Analytical BS'), 
title('Variance Gamma Model')
legend
axis([min(strikes) max(strikes) 0.0 1800])
subplot(2,4,7)
plot(K,(C_COS_vg-C_BS')), grid on, hold on;
title("Diff. between Variance Gamma" + newline + "and Analytical Black-Scholes")
axis([min(strikes) max(strikes) -5.0 70.0])

% CGMY
subplot(2,4,4)
plot(K,C_COS_cgmy,'r', 'DisplayName', 'CGMY'), grid on, hold on;
plot(K,C_BS,'--k', 'DisplayName', 'Analytical BS'), 
title('CGMY Model')
legend
axis([min(strikes) max(strikes) 0.0 1800])
subplot(2,4,8)
plot(K,(C_COS_cgmy-C_BS')), grid on, hold off;
title("Diff. between CGMY" + newline + "and Analytical Black-Scholes")
axis([min(strikes) max(strikes) -5.0 70.0])



%% ========== PRIVATE FUNCTIONS ==========
function MakePdfPlot2(cf,a,b,k,mu,u_0,T,N,plotpos,name)
%  Based on Peter.Gruber@unisg.ch, February 2007
bma = b-a;
x = linspace(a,b,N);
true=normpdf(x,mu, sqrt(u_0 * T));

Fk = 2/bma * real( cf.*exp(-1i*k*a*pi/bma)   );
Fk(1)=0.5*Fk(1);

for l=1:length(x)
    pdf(l)=sum(Fk.*cos (k*pi*(x(l)-a)/bma) );
end

subplot(2,4,plotpos)
hold on
plot(x,pdf, 'DisplayName', name), grid on;
title(name)
plot(x, true, 'DisplayName', 'Normal Pdf')
legend
axis([-sqrt(u_0*T)*12 sqrt(u_0*T)*12 0 4])
hold off

subplot(2,4,plotpos+4)
plot(x,pdf(:)-true(:)), grid on;
title(strcat(name, " - Normal Pdf"))
axis([-sqrt(u_0*T)*12 sqrt(u_0*T)*12 -1 1.2])

end
