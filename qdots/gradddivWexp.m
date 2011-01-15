function gf = gradddivWexp(Wxk_r, varargin)
% function gf = gradddivW(Wxk, varargin)
% gf is 1 x x*k row vector
% Wxk_r = reshape(Wxk,1,x*k) -> row vector
% Vxt = varargin{1};  %data
% Hkt = varargin{2};  %H matrix
% Wxk_fix = varargin{3}; %fixed part of the Wxk matrix (e.g. background) ->rows
% Hkt_fix = varargin{4}; %fixed part (lines) of the H matrix (e.g. background)

alphaW=1; %for now....

Vxt = varargin{1};  %data
Hkt_tmp = varargin{2};  %H matrix
Wxk_fix = varargin{3}; %fixed part of the Wxk matrix (e.g. background) ->rows
Hkt_fix = varargin{4}; %fixed part (lines) of the H matrix (e.g. background)
peval = varargin{5}; %parameters

if ~isfield(peval, 'w_lambda') peval.w_lambda=0; end

x=size(Vxt,1);
k=length(peval.w_dovec);

Wxk_tmp = exp(reshape(Wxk_r,x,k));

Wxk = zeros(peval.numpix, peval.ncomp+1);
Hkt = zeros(peval.ncomp+1, peval.nt);

Wxk(:,peval.w_dovec)=Wxk_tmp;
Hkt(peval.h_dovec,:)=Hkt_tmp;

Wxk(:,peval.w_fixvec)=Wxk_fix;
Hkt(peval.h_fixvec,:)=Hkt_fix;

sHxk = repmat((sum(Hkt,2))',size(Vxt,1),1);
gfxk = (sHxk - (Vxt./(Wxk*Hkt))*Hkt')*alphaW.*Wxk; %d/dw(d-divergence)
gf = reshape(gfxk(:,1:k),1,x*k);