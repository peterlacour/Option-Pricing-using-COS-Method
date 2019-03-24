clear, clc;

% Default Parameters Fang and Oosterlee (2008)
lambda      =  1.5768;                          % Mean Reversion Speed
u_bar       =  0.0398;                          % Long-run Variance
u_0         =  0.0175;                          % Initial Variance
eta         =  0.5751;                          % Volatility of the Volatility Parameter
rho         = -0.5711;                          % Correlation between Wiener Processes


% Option / Underlying / Market Parameters
T    = 1.00;                                    % Time to Maturity
r    = 0.0;                                       % Risk-free Rate
q    = 0;                                       % Dividend Yield
mu   = r - q;                                   % Price Drift Rate
S0   = 100;

%% SETUP
Npoints = 5000;
x = linspace(-2,2,Npoints);
N = 5000;                         % default value 25, try 100,500,5000
k = 0:N-1;

%{
% Heston Plot
% Heston truncation bounds
[a, b] = cos_truncation_range_v2(mu, lambda, u_bar, u_0, eta, rho, T);
phi_hes =@(t) heston_char_function_v1(mu, lambda, u_bar, u_0, eta, rho, t, T);
x = linspace(a,b,Npoints);
true=normpdf(x,mu, sqrt(u_bar) );
figure
MakePdfPlot2(phi_hes,a,b,x,true,N,1, 'Heston Pdf');
%}

x = linspace(-2,2,Npoints);
% Variance Gamma Plot
theta = -0.1436;
[c1, c2, c4, ~] = variance_gamma_cumulants_v2( u_0, T, theta, mu );
[a, b]    = cos_truncation_range_v2(c1,c2,c4,10);
phi_vg = vg_char_fn(u_0, theta, a, b, k, T);
x = linspace(a,b,Npoints);
true=normpdf(x,mu, sqrt(u_0) );
MakePdfPlot2(phi_vg,a,b,x,true,N,2, 'Variance Gamma Pdf');

%{
% CGMY Plot
% CGMY - Parameters from CGMY ( 2003 )
C = 24.79;
G = 94.45;
Y = 0.2495;
M = 95.79;
[a, b] = cos_cgmy_truncation_range_v1(u_bar, T, mu, C, G, M, Y);
x = linspace(a,b,Npoints);
phi_cgmy =@(t) cgmy_char_fn(t, u_bar, T, mu, C, G, Y, M);
MakePdfPlot2(phi_cgmy,a,b,x,true,N,3, 'CGMY Pdf');


% FMLS Plot
alpha = 1.1;
phi_fmls =@(t) fmls_char_fn(alpha, t, eta, T, mu);
MakePdfPlot2(phi_fmls,-2,2,x,true,N,4, 'FMLS Pdf');

%}

% BS Plot
[c1, c2, c4, ~] = bs_cumulants_v1(u_0, mu, T );
[a, b]    = cos_truncation_range_v2(c1,c2,c4,10);
x = linspace(a,b,Npoints);
phi_bs = bs_char_fn_v1(mu, u_0, a, b, k, T);
MakePdfPlot2(phi_bs,a,b,x,true,N,3, 'BS Pdf');



%% ========== PRIVATE FUNCTIONS ==========
function MakePdfPlot2(cf,a,b,x,true,N,plotpos,name)
bma = b-a;
k = 0:N-1;
Fk = 2/bma * real( cf .*exp(-1i*k*a*pi/bma)   );
Fk(1)=0.5*Fk(1);
for l=1:length(x)
    pdf(l)=sum(Fk.*cos (k*pi*(x(l)-a)/bma) );
end

subplot(2,3,plotpos)
hold on
plot(x,pdf, 'DisplayName', name);
title(name)
%a=axis;
%axis([min(x) max(x) a(3:4)])
plot(x, true, 'DisplayName', 'Normal Pdf')
legend
hold off

subplot(2,3,plotpos+3)
plot(x,pdf(:)-true(:))
title(strcat(name, " - Normal Pdf"))
%a=axis;
%axis([min(x) max(x) a(3:4)])
%end

end