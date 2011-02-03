function f = ddivW(Wxk_r, varargin)
% function f = ddivW(Wxk_r, varargin)
% Wxk_r = reshape(Wxk,1,x*k) -> row vector
% Vxt = varargin{1};  %data
% Hkt = varargin{2};  %H matrix

Vxt = varargin{1};  %data
Hkt = varargin{2};  %H matrix
x=size(Vxt,1);
k=size(Hkt,1);
Wxk=reshape(Wxk_r,x,k);

fxt = Vxt.*log(Vxt./(Wxk*Hkt))-Vxt+Wxk*Hkt; %d-divergence 
f = sum(fxt(:));
