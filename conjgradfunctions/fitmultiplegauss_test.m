function f=fitmultiplegauss_test(x, varargin)
% f=loglikGaP(x, varargin)
% complete log likellihood of the GaP model 
% Vxt = varargin{1};      %data
% sigpsf = varargin{2};  %std deviation of the PSF gaussian approx
% alpha = varargin{3}; %parameters of the Gamma prior on the blinking
% beta = varargin{4}; %parameters of the Gamma prior on the blinking
% peval = varargin{5}; %parameters
% x(1:end-2*peval.ncomp) is Hkt


Vxyt = varargin{1};      %data
sigpsf = varargin{2};  %std deviation of the PSF gaussian approx
% alpha = varargin{3}; %parameters of the Gamma prior on the blinking
% beta = varargin{4}; %parameters of the Gamma prior on the blinking
peval = varargin{3}; %parameters

Vxt=reshape(Vxyt,peval.nx*peval.ny,1);
sigpsf_vec=repmat(sigpsf,peval.ncomp,1); %all psfs same sigma
a_vec=1./(sigpsf_vec.^2*2*pi); % all normalised to 1
cxy_vec=[x(1:peval.ncomp)'+1, x(peval.ncomp+1:end)'+1]; %different notation of the dip_image/new version

Wxkpix=gauss2dmultislice([peval.nx, peval.ny, peval.ncomp], cxy_vec, sigpsf_vec, a_vec);
Wxkpix=normalize(Wxkpix); %normalize PSFs to 1
Wxk=reshape(Wxkpix,peval.nx*peval.ny, peval.ncomp);

factor=(sum(Vxt(:))/sum(Wxk(:)))/peval.ncomp;
Hkt=factor*ones(peval.ncomp,1); %just sum the components


P=Wxk*Hkt+peval.bg; %current approximation

%Poisson contribution
t1=Vxt.*log(P) - P;
%Gamma contribution
% t2=(alpha-1)*log(Hkt)-1/beta*Hkt-alpha*log(beta)-log(gamma(alpha));

f=sum(t1(:));%+sum(t2(:)));
f=-f; %conjugate gradient is mimimizing!
end