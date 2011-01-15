function gf=gradloglikGaP(x, varargin)
% f=loglikGaP(x, varargin)
% complete log likellihood of the GaP model 
% Vxt = varargin{1};      %data
% sigpsf = varargin{2};  %std deviation of the PSF gaussian approx
% alpha = varargin{3}; %parameters of the Gamma prior on the blinking
% beta = varargin{4}; %parameters of the Gamma prior on the blinking
% peval = varargin{5}; %parameters
% x(1:end-2*peval.ncomp) is Hkt


Vxt = varargin{1};      %data
sigpsf = varargin{2};  %std deviation of the PSF gaussian approx
alpha = varargin{3}; %parameters of the Gamma prior on the blinking
beta = varargin{4}; %parameters of the Gamma prior on the blinking
peval = varargin{5}; %parameters

[Hkt, cx, cy, Wxk] = loglikGaPreadparam(x,varargin);


[Wxkbg,Hktbg]=addbg(Wxk, Hkt, peval.bg);
Hkt=max(Hkt,eps); % adjust small values to avoid undeflow
P=Wxkbg*Hktbg; %current approximation

%linear grasdient shifted by cx
xxvc = lineargrad([peval.nx, peval.ny, peval.ncomp], cx, 'xx'); 
yyvc = lineargrad([peval.nx, peval.ny, peval.ncomp], cy, 'yy');

% dW/dcx:
Wxtcx=1/sigpsf^2*xxvc.*Wxk; 
% dW/dcy:
Wxtcy=1/sigpsf^2*yyvc.*Wxk;

% d(log(L))/dHkt:
% gfHkt=(alpha-1)*1./Hkt - 1/beta*Hkt + Wxk'*(Vxt./P)-1;
gfHkt= Wxk'*(Vxt./P)-1; %without background (->not Wxkgb) and d(log(L)/dcx)
% d(log(L))/dcx:
gfcx=diag(Wxtcx'*(Vxt./P-ones(peval.nx*peval.ny, peval.nt))*Hkt');
% d(log(L))/dcy:
gfcy=-diag(Wxtcy'*(Vxt./P-ones(peval.nx*peval.ny, peval.nt))*Hkt'); 

% gf = [reshape(gfHkt,1,peval.nt*peval.ncomp), gfcx', gfcy'];
gf = [reshape(gfHkt,1,peval.nt*peval.ncomp)];
% gf = [gfcx', gfcy'];
gf=-gf; %conjugate gradient is minimizing!
end
function Wnorm=normalizePSF(W)
sw=size(W);
Wr=reshape(W, sw(1)*sw(2),sw(3));
q=squeeze(sum(Wr,1));
Wrnorm=Wr./repmat(q,sw(1)*sw(2),1);
Wnorm=reshape(Wrnorm,sw(1), sw(2), sw(3));
end
function xxvc = lineargrad(sizevec, cx, dir)
switch dir
    case 'xx'
%         xxp=double(xx(sizevec, 'true')); %linear function - pixels
        xxp=double(xx(sizevec, 'corner')); %linear function - pixels
    case 'yy'
%         xxp=double(yy(sizevec, 'true')); %linear function - pixels
        xxp=double(yy(sizevec, 'corner')); %linear function - pixels
    otherwise 
        error('Wrong dir')        
end
    
xxv=reshape(xxp,sizevec(1)*sizevec(2),sizevec(3)); %linear function - vector
xxvc=xxv-repmat(cx,sizevec(1)*sizevec(2),1);
% xxp=double(xx([peval.nx, peval.ny, peval.ncomp], 'true')); %linear function - pixels
%yyp=double(yy([peval.nx, peval.ny, peval.ncomp], 'true'));
% xxv=reshape(xxp,peval.nx*peval.ny,peval.ncomp); %linear function - vector
%yyv=reshape(yyp,peval.nx*peval.ny,peval.ncomp);
% xxvc=xxv-repmat(cx,peval.ncomp,1);
%yyvc=yyv-repmat(cy,peval.ncomp,1);
end