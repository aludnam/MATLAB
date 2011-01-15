function f = ddivH(Hkt_r, varargin)
% function f = ddivH(Hkt_r, varargin)
% Hkt_r = reshape(Hkt,1,k*t) -> row vector
% Vxt = varargin{1};  %data
% Wxk = varargin{2};  %W matrix

Vxt = varargin{1};  %data
Wxk = varargin{2};  %W matrix
t=size(Vxt,2);
k=size(Wxk,2);
Hkt=reshape(Hkt_r,k,t);

fxt = (Vxt.*log(Vxt./(Wxk*Hkt))-Vxt+Wxk*Hkt); %d-divergence 
f = sum(fxt(:));
