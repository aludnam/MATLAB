function nout=gammalogapprox(nin,method)
% Computer approximation of the log-gamma (log(gamma(n)))
% from: http://mathworld.wolfram.com/StirlingsApproximation.html
% method:   'stirling' Striling approximation
%           'Gosper'  (default) Gosper approximation (better?)

if ~exist('method', 'var')
    method = 'Gosper'; % Default
end

nout = factorialapprox(nin-1, method);

