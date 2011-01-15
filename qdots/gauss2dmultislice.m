function G=gauss2dmultislice(sizevec, m, s, a)
% G=gauss2dmultislice(sizevec, m, s, a)
lm=length(m);
if length(s)<lm
    s=repmat(s,lm); %all the sama sigma
end
if length(a)<lm
    a=repmat(a,lm); %all the sama amplitude
end

G=zeros(sizevec);
for ii=1:sizevec(3)    
    G(:,:,ii) = gauss2d(sizevec(1:2), m(ii,:), s(ii), a(ii));
end