function gf=gradfitmultiplegauss_test(x, varargin)
% f=loglikGaP(x, varargin)
% complete log likellihood of the GaP model 
% Vxt = varargin{1};      %data
% sigpsf = varargin{2};  %std deviation of the PSF gaussian approx
% alpha = varargin{3}; %parameters of the Gamma prior on the blinking
% beta = varargin{4}; %parameters of the Gamma prior on the blinking
% peval = varargin{5}; %parameters
% x(1:end-2*peval.ncomp) is Hkt

Vxyt = varargin{1};      %data
sigpsf = varargin{2};  %std deviation of the PSF gaussian approx
% alpha = varargin{3}; %parameters of the Gamma prior on the blinking
% beta = varargin{4}; %parameters of the Gamma prior on the blinking
peval = varargin{3}; %parameters

Vxt=reshape(Vxyt,peval.nx*peval.ny,1);
sigpsf_vec=repmat(sigpsf,peval.ncomp,1); %all psfs same sigma
a_vec=1./(sigpsf_vec.^2*2*pi); % all normalised to 1
cxy_vec=[x(1:peval.ncomp)'+1, x(peval.ncomp+1:end)'+1]; %different notation of the dip_image/new version

Wxkpix=gauss2dmultislice([peval.nx, peval.ny, peval.ncomp], cxy_vec, sigpsf_vec, a_vec);
Wxkpix=normalize(Wxkpix); %normalize PSFs to 1
Wxk=reshape(Wxkpix,peval.nx*peval.ny, peval.ncomp);

factor=(sum(Vxt(:))/sum(Wxk(:)))/peval.ncomp;
Hkt=factor*ones(peval.ncomp,1); %just sum the components


P=Wxk*Hkt+peval.bg; %current approximation

%linear grasdient shifted by cx
xxvc = lineargrad([peval.nx, peval.ny, peval.ncomp], cxy_vec(:,1)', 'xx'); 
yyvc = lineargrad([peval.nx, peval.ny, peval.ncomp], cxy_vec(:,2)', 'yy');

% dW/dcx:
Wxtcx=1/sigpsf^2*xxvc.*Wxk; 
% dW/dcy:
Wxtcy=1/sigpsf^2*yyvc.*Wxk;

% d(log(L))/dHkt:
% gfHkt=(alpha-1)*1./Hkt - 1/beta*Hkt + Wxk'*(Vxt./P)-1;
gfHkt= Hkt.*(Wxk'*(Vxt./P)-1); %without background (->not Wxkgb) and d(log(L)/dcx)
% d(log(L))/dcx:
gfcx=diag(Wxtcx'*(Vxt./P-ones(peval.nx*peval.ny, 1))*Hkt');
% d(log(L))/dcy:
gfcy=diag(Wxtcy'*(Vxt./P-ones(peval.nx*peval.ny, 1))*Hkt'); 

% gf = [reshape(gfHkt,1,peval.nt*peval.ncomp), gfcx', gfcy'];
% gf = [reshape(gfHkt,1,peval.nt*peval.ncomp)];
gf = [gfcx', gfcy'];
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