tau=c;
s=.1;
a=(1./tau)*s;
N=1/s;


% m=N/(exp(a*N)-1)-exp(a)/(exp(a)-1);
% this is the mean: sum_n(a*n*exp(-a*n))
m=-s*exp(-a*N)*(N*(exp(a)-1)-exp(a)*(exp(a*N)-1))/(exp(a)-1)^2

%this is normalization constant
q=(exp(a)-exp(-a*N))/(exp(a)-1)