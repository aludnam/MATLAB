function d = ddivergence(A,B)
% d = ddivergence(A,B);
dm = bsxfun(@times, A, log (bsxfun(@rdivide,A,B))) - A + B;
d = sum(dm(:));
% d = sum(sum(A .* log (A./B))); %KLdiv
end