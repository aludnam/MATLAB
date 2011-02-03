function f = ddivWexp(Wxk_r, varargin)
% function f = ddivW(Wxk_r, varargin)
% Wxk_r = reshape(Wxk,1,x*k) -> row vector
% Vxt = varargin{1};        %data
% Hkt = varargin{2};        %H matrix
% Wxk_fix=varargin{3};      %fixed part of the Wxk matrix (e.g. background)
% Hkt_fix = varargin{4};    %fixed part (lines) of the H matrix (e.g. background)

Vxt = varargin{1};      %data
Hkt_tmp = varargin{2};  %H matrix
Wxk_fix = varargin{3};  %fixed part of the Wxk matrix (e.g. background) ->rows
Hkt_fix = varargin{4};  %fixed part (lines) of the H matrix (e.g. background)
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

%%%%
% sumw = sum(Wxk,1);
% Wxk = Wxk./repmat(sumw,size(Wxk,1),1);


% fxt = Vxt.*log(Vxt./(Wxk*Hkt))-Vxt+Wxk*Hkt; %d-divergence 
% f = sum(fxt(:));
% peval.w_lambda*Wxk.^2
f = ddivergence(Vxt,Wxk*Hkt);