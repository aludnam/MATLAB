function nout=gammalogapprox(nin,method)
% Compute approximation of the log-gamma (log(gamma(n)))
% from: http://mathworld.wolfram.com/StirlingsApproximation.html
% method:   'Stirling' Striling approximation
%           'Gosper'  (default) Gosper approximation (better?)

if ~exist('method', 'var')
    method = 'Gosper'; % Default
end

if nin-1==0
    nout=1;
else
    nout = factorialapprox(nin-1, method);
end

