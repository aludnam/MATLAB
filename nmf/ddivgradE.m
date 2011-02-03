function gf = ddivgradE(Ekt, varargin)
% function f = ddivE(Ekt, varargin)
% Vxt = varargin{1};  %data
% Wxk = varargin{2};  %estimated H from NMF update
% alpha = varargin{3}; %rate

Vxt = varargin{1};  %data
Wxk = varargin{2};  %estimated H from NMF update
alpha = varargin{3};  %rate
nt=size(Vxt,2);
nk=size(Wxk,2);

Hkt=exp(alpha*reshape(Ekt,nt,nk));

x1=repmat(sum(Wxk,1)',1,nt);
gfmat=alpha*Hkt.*(Wxk'*(Vxt./(Wxk*Hkt))-x1);

gf=reshape(gfmat',nk*nt,1); %nk*nt x 1 vector


