function gf=gradloglikGaPExp(x, varargin)
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

[Hkt_r, cx, cy, Wxk] = loglikGaPreadparam(x,varargin);

Hkt=exp(Hkt_r); %nonnegativity constrains
[Wxkbg,Hktbg]=addbg(Wxk, Hkt, peval.bg);

P=Wxkbg*Hktbg; %current approximation

%linear grasdient shifted by cx
xxvc = lineargrad([peval.nx, peval.ny, peval.ncomp], cx, 'xx'); 
yyvc = lineargrad([peval.nx, peval.ny, peval.ncomp], cy, 'yy');

% dW/dcx:
Wxtcx=1/sigpsf^2*xxvc.*Wxk; 
% dW/dcy:
Wxtcy=1/sigpsf^2*yyvc.*Wxk;

% d(log(L))/dHkt:
% gfHkt=(alpha-1)*1./Hkt - 1/beta + Wxk'*(Vxt./P)-1;
gfHkt= Hkt.*((Wxk'*(Vxt./P)-1)-(alpha-1)./Hkt - ones(size(Hkt))*1/beta); %without background (->not Wxkgb) and d(log(L)/dcx)
% d(log(L))/dcx:
gfcx=diag(Wxtcx'*(Vxt./P-ones(peval.nx*peval.ny, peval.nt))*Hkt');
% d(log(L))/dcy:
gfcy=diag(Wxtcy'*(Vxt./P-ones(peval.nx*peval.ny, peval.nt))*Hkt'); 

gf = [reshape(gfHkt,1,peval.nt*peval.ncomp), gfcx', gfcy'];
% gf = [reshape(gfHkt,1,peval.nt*peval.ncomp)];
% gf = [gfcx', gfcy'];
gf=-gf; %conjugate gradient is minimizing!
end

function xxvc = lineargrad(sizevec, cx, dir)
switch dir
    case 'xx'
        xxp=double(xx(sizevec, 'corner')); %linear function - pixels
    case 'yy'
        xxp=double(yy(sizevec, 'corner')); %linear function - pixels
    otherwise 
        error('Wrong dir')        
end    
xxv=reshape(xxp,sizevec(1)*sizevec(2),sizevec(3)); %linear function - vector
xxvc=xxv-repmat(cx,sizevec(1)*sizevec(2),1);
end