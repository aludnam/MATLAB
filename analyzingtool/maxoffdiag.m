function v=maxoffdiag(A)
% v=maxoffdiag(A)
% Maximum non-diagonal element of abs(A). If ndims(A)=3 then v is a vector
% with max non-diag element of each slice (along 3ddimension) of A.

ns = size(A,3);
v=zeros(1,ns);
for ii=1:ns
    An = (A(:,:,ii) - diag(diag(A(:,:,ii))));
    v(ii) = max(abs(An(:)));
end

