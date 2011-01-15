function gf = gradddivW(Wxk_r, varargin)
% function gf = gradddivW(Wxk, varargin)
% gf is 1 x x*k row vector
% Wxk_r = reshape(Wxk,1,x*k) -> row vector
% Vxt = varargin{1};  %data
% Hkt = varargin{2};  %H matrix

Vxt = varargin{1};  %data
Hkt = varargin{2};  %H matrix
x=size(Vxt,1);
k=size(Hkt,1);
Wxk=reshape(Wxk_r,x,k);

sHxk = repmat((sum(Hkt,2))',size(Vxt,1),1);
gfxk = sHxk - (Vxt./(Wxk*Hkt))*Hkt'; %d/dw(d-divergence)
gf = reshape(gfxk,1,x*k);