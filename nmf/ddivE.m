function f = ddivE(Ekt, varargin)
% function f = ddivE(Ekt, varargin)
% Vxt = varargin{1};  %data
% Wxk = varargin{2};  %estimated H from NMF update
% alpha = varargin{3}; %rate

Vxt = varargin{1};  %data
Wxk = varargin{2};  %estimated values of W
alpha = varargin{3};

Hkt=exp(alpha*Ekt);
f = ddivergence(Vxt,Wxk*Hkt);

