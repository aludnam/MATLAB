function gf = gradspars(Hkt_r, varargin)
% function gf = gradddivH(Hkt_r, varargin)
% gf is 1 x k*t row vector
% Hkt_r = reshape(Hkt,1,k*t) -> row vector
% Vxt = varargin{1};  %data
% Wxk = varargin{2};  %W matrix
% Wxk_fix = varargin{3}; %fixed part of the Wxk matrix (e.g. background) ->rows
% Hkt_fix = varargin{4}; %fixed part (lines) of the H matrix (e.g. background)

alphaH=1; %for now....
Vxt = varargin{1};  %data
Wxk_tmp = varargin{2};  %W matrix
Wxk_fix = varargin{3}; %fixed part of the Wxk matrix (e.g. background) ->rows
Hkt_fix = varargin{4}; %fixed part (lines) of the H matrix (e.g. background)
peval = varargin{5}; %parameters

if ~isfield(peval, 'w_lambda') peval.w_lambda=0; end

t=size(Vxt,2);
k=length(peval.h_dovec);

Hkt_tmp = exp(reshape(Hkt_r,k,t));

Wxk = zeros(peval.numpix, peval.ncomp+1);
Hkt = zeros(peval.ncomp+1, peval.nt);

Wxk(:,peval.w_dovec)=Wxk_tmp;
Hkt(peval.h_dovec,:)=Hkt_tmp;

Wxk(:,peval.w_fixvec)=Wxk_fix;
Hkt(peval.h_fixvec,:)=Hkt_fix;

deltasum=sum(sum(Wxk_tmp,1))-k;
if and(~isempty(Wxk_tmp), abs(deltasum)>10^-6)
    error('Wxk is not correctly normalized! (sum(Wxk_tmp,1)<>1)\n sum(Wxk_tmp,1)=%f',deltasum)
end


% ap = peval.alphapen;
% sumH_t = sum(Hkt(peval.h_dovec,:),2);
sumH_t = sum(Hkt,2);
sumH_t_sq = sum(sumH_t.^2);
sumH = sum(sumH_t);
nh = length(sumH_t); 

sparsness = (sqrt(nh)-sumH/(sqrt(sumH_t_sq)))/(sqrt(nh)-1);
% gradspars=(1/(sqrt(nh)-1))*(sumH/(sumH_t_sq)^1.5*sumH_t - 1/(sqrt(sumH_t_sq)));
gradspars= (1/(sqrt(nh)-1)*(sumH_t*(sumH)-sumH_t_sq)/(sumH_t_sq^1.5));
gradspars_kt = repmat(gradspars, 1, t);

% fkt = ddivergence(Vxt,Wxk*Hkt);
% gfkt = (1-Wxk'*(Vxt./(Wxk*Hkt)))*alphaH.*Hkt; %d/dh(d-divergence)
% gfkt = (1-Wxk'*(Vxt./(Wxk*Hkt))).*Hkt; %d/dh(d-divergence)
%(fg)'=(f'g-fg')/g^2;
% gfkt_all=(gfkt.*sparsness-fkt.*gradspars_kt)./(sparsness).^2;
% gfkt_all = gradspars_kt.*Hkt(peval.h_dovec, :); %because of exp...
gfkt_all = gradspars_kt.*Hkt; %because of exp...
% one is tehre because Wxt is normalized: sum(Wxt,1)=1
gf=reshape(gfkt_all(peval.h_dovec,:),1,k*t); %making row vector