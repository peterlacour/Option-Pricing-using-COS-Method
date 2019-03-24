function  [c,p,d1,d2] = black_scholes_price(S,X,r,T,sigma,q)
%Calculates Black-Scholes european option prices.
%
%  Usage:      [c,p,d1,d2] = blackS( S,X,r,T,sigma,[q] );
%
%  Inputs:     S      scalar or nx1 vector, possible current stock prices
%              X      scalar or nx1 vector, strike price
%              r      scalar, riskfree interest rate (continuously compounded)
%              T      scalar, time to expiry of option
%              sigma  scalar or nx1 vector, std in stock price evolution
%              [q]    scalar, dividend yield (continuously compounded), optional
%
%  Output:     c      nx1 vector, call option prices
%              p      nx1 vector, put option prices
%              d1     nx1 vector
%              d2     nx1 vector
%
%  Peter.Gruber@unisg.ch, February 2007
%  Based on code by Paul.Soderlind@unisg.ch
if nargin==6     % if dividend is specified, correct for it
    S = S * exp(-q*T);
end
d1 = ( log(S./X) + (r + 1/2*sigma.^2)*T ) ./ (sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);
c  = S.*stdnCdf(d1) - X.*exp(-r*T).*stdnCdf(d2);
p  = c + X.*exp(-r*T) - S;                  %put-call parity
end

function cdf = stdnCdf(a)
cdf = 0.5 + 0.5*erf(a/sqrt(2));
end

