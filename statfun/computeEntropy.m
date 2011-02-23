function p=computeEntropy(in)
% p=computeEntropy(in)
% Computer entropy of the (nomarlized) vector

innorm = max(in./sum(in), eps); %to avoid 0 
p=-sum(innorm.*log(innorm));