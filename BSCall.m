function bsc = BSCall(s,k,sigma,tau,r)

% BSCALL computes the price of a Black-Scholes European call option.
%
% Inputs: k = strike price
%         s = price of underlying asset
%         sigma = annual volatility
%         tau = time to maturity (as a fraction of a year)
%         r = interest rate
%

echo off;
stau=sqrt(tau);
d1=(log(s./(k.*exp(-r.*tau))))./(sigma.*stau)+0.5*sigma.*stau;
d2=d1-sigma.*stau;
bsc=s.*normcdf(d1,0,1)-k.*exp(-r.*tau).*normcdf(d2,0,1);
echo on;

end