function d = loglikelihoodPoisson(A,B)
% d = loglikelihoodPoisson(A,B)
% For NMF: A=data, B=model
    dm = A .* log (B) - B - (A.*log(A)-A); %stirling approximation of log(A!)
    d = sum(dm(:));
end