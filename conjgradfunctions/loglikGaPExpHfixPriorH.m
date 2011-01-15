function f=loglikGaPExpHfixPriorH(x, varargin)
% f=loglikGaP(x, varargin)
% complete log likellihood of the GaP model 
% Vxt = varargin{1};      %data
% sigpsf = varargin{2};  %std deviation of the PSF gaussian approx
% alpha = varargin{3}; %parameters of the Gamma prior on the blinking
% beta = varargin{4}; %parameters of the Gamma prior on the blinking
% peval = varargin{5}; %parameters
% x(1:end-2*peval.ncomp) is Hkt


Vxt = varargin{1};      %data
sigpsf = varargin{2};  %std deviation of the PSF gaussian approx
alpha = varargin{3}; %parameters of the Gamma prior on the blinking
beta = varargin{4}; %parameters of the Gamma prior on the blinking
peval = varargin{5}; %parameters

[Hkt_r, cx, cy, Wxk] = loglikGaPreadparamHfix(x,varargin);
Hkt=exp(Hkt_r); %nonnegativity constrains

[Wxkbg,Hktbg]=addbg(Wxk, Hkt, peval.bg);
P=Wxkbg*Hktbg; %current approximation

%Poisson contribution
t1=Vxt.*log(P) - P;
%Gamma contribution
t2=(alpha-1)*log(Hkt)-1/beta*Hkt-alpha*log(beta)-log(gamma(alpha));

f=sum(t1(:))+sum(t2(:));
f=-f; %conjugate gradient is mimimizing!
end