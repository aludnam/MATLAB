function gf = gradddivHexpS51(Hkt_r, varargin)
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


t=size(Vxt,2);
k=size(Wxk_tmp,2);

Wxk = [Wxk_tmp, Wxk_fix];
Hkt=[exp(reshape(Hkt_r,k,t)); Hkt_fix];

deltasum=sum(sum(Wxk_tmp,1))-k;
if abs(deltasum)>1e-6
    error('Wxk is not correctly normalized! (sum(Wxk_tmp,1)<>1)\n sum(Wxk_tmp,1)=%f',deltasum)
end

gfkt = (1-Wxk'*(Vxt./(Wxk*Hkt)))*alphaH.*Hkt; %d/dh(d-divergence)
% one is tehre because Wxt is normalized: sum(Wxt,1)=1
gf=reshape(gfkt(1:k,:),1,k*t); %making row vector