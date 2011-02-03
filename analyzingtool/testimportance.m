function testscore = testimportance(V,W,H) 

k = size(W,2);
for ii=1:k
    Wn = (removerows(W',ii))'; %remove columns
    Hn = removerows(H,ii); 
    testscore(ii)=ddivergence(V,Wn*Hn);
end
    