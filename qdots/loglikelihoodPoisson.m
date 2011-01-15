function d = loglikelihoodPoisson(A,B)
% d = ddivergence(A,B);
    dm = A .* log (B) - B - (A.*log(A)-A); %stirling approximation of log(A!)
    d = sum(dm(:));
end