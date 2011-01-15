function d = Ddiv(A,B)
    d = sum(sum(A .* log (A./B) - A + B));
end