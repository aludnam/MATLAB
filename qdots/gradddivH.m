function gf = gradddivH(Hkt_r, varargin)
% function gf = gradddivH(Hkt_r, varargin)
% gf is 1 x k*t row vector
% Hkt_r = reshape(Hkt,1,k*t) -> row vector
% Vxt = varargin{1};  %data
% Wxk = varargin{2};  %W matrix

Vxt = varargin{1};  %data
Wxk = varargin{2};  %W matrix
t=size(Vxt,2);
k=size(Wxk,2);
Hkt=reshape(Hkt_r,k,t);

if ~(sum(sum(Wxk,1))==k)
    error('Wxk is not correctly normalized! (sum(Wxk,1)<>1)')
end

gfkt = 1-Wxk'*(Vxt./(Wxk*Hkt)); %d/dh(d-divergence)
% one is tehre because Wxt is normalized: sum(Wxt,1)=1
gf=reshape(gfkt,1,t*k); %making row vector