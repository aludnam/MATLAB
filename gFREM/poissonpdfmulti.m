function A=poissonpdfmulti(x,lambda_vec)

A=zeros(length(x), size(lambda_vec,1));
for ii=1:size(lambda_vec,1)
    A(:,ii)=poisspdf(x,lambda_vec(ii));
end