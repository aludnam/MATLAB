function [Hkt, cx, cy, Wxk] = loglikGaPreadparamWfix(x,qqq)
% [Hkt, cx, cy, Wxk] = loglikGaPreadparam(varargin);

sigpsf = qqq{2};  %std deviation of the PSF gaussian approx
peval =  qqq{5}; %parameters
% Hkt_linear=x(1:end-peval.ncomp*2); % intensities
% cx=x(end-peval.ncomp*2+1:end-peval.ncomp); %x-coordinates of the centers
% cy=x(end-peval.ncomp+1:end); % y-coordinates of the centers
% Hkt_linear=qqq{6}; % intensities
Hkt_linear=x; % intensities
cx= qqq{6};
cy= qqq{7};

sigpsf_vec=repmat(sigpsf,peval.ncomp,1); %all psfs same sigma

% this is anoying as it depends when the simulation was performed/... grr
cxy_vec=[cx'+1, cy'+1]; %different notation of the dip_image/new version
% cxy_vec=[cx'+0.5, cy'+0.5]; %different notation of the dip_image/old version
a_vec=1./(sigpsf_vec.^2*2*pi); % all normalised to 1

Hkt=reshape(Hkt_linear, peval.ncomp, peval.nt);
% generate PSFs from given parameters:
Wxkpix=gauss2dmultislice([peval.nx, peval.ny, peval.ncomp], cxy_vec, sigpsf_vec, a_vec);
Wxkpix=normalizePSF(Wxkpix); %normalize PSFs to 1
Wxk=reshape(Wxkpix,peval.nx*peval.ny, peval.ncomp);
end

function Wnorm=normalizePSF(W)
sw=size(W);
if numel(sw)<3 %if ncomp=1... 
    sw(3) = 1;
end
Wr=reshape(W, sw(1)*sw(2),sw(3));
q=squeeze(sum(Wr,1));
Wrnorm=Wr./repmat(q,sw(1)*sw(2),1);
Wnorm=reshape(Wrnorm,sw(1), sw(2), sw(3));
end
