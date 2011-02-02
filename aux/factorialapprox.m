function nout = factorialapprox(nin, method)
% Compute approximation of the log-factorial (log(nin!))
% from: http://mathworld.wolfram.com/StirlingsApproximation.html
% method:   'Stirling' Striling approximation
%           'Gosper'  (default) Gosper approximation (better?)

if ~exist('method', 'var')
    method = 'Gosper'; % Default
end

if nin==0
    nout=1;
else
    switch lower(method)
        case 'stirling'
            nout=(nin+.5).*log(nin) - nin + .5*log(2*pi);
        case 'gosper'
            nout=nin.*log(nin)-nin+0.5*log((2*nin+1/3)*pi);
    end
end


