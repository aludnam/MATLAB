function [varargout, peval]=updates_nmfclassic(peval,varargin)
% [varargout, peval]=updates_nmfclassic(peval,varargin)
% Non-negative matrix factorisation updates minimising KL divergence
% KL(V|WH) 
% references : D.D. Lee and H.S. Seung. Algorithms for non-negative matrix
% factorization. Advances in neural information processing systems, 13, 2001.
%
% v    : N x T data matrix
% w    : N x K NMF factor
% h    : K x T NMF factor
% K is the rank of factorisation. 
%
% v = varargin{1};
% w = varargin{2};  % initial values  
% h = varargin{3};  % initial values
% peval.h_fixvec & peval.w_fixvec:
% fixvec - vector which component should be fixed: eg [2 3] will
% fix second and third component while varying the first...

mfprintf(peval.fid,'Classic NMF iterations.\n')
peval=setDefaultValuesPeval(peval);
v = varargin{1};
w = varargin{2};
h = varargin{3};
[peval.dovec_w, peval.dovec_h]=setDoVec(peval);

checkv(v)           %check values of v
checkddivfreq=50;   % How often ot check ddiv
dall=zeros(1, ceil(peval.maxiter/checkddivfreq));
indexd=1;
for ii=2:peval.maxiter    
    
    y1=w'*(bsxfun(@rdivide,v,w*h));    
    sw = sum(w,1)';
    h(peval.dovec_h,:)= bsxfun(@rdivide,h(peval.dovec_h,:).*(y1(peval.dovec_h,:)),sw(peval.dovec_h));
    h=max(h,eps); % adjust small values to avoid undeflow
    
    y2=(bsxfun(@rdivide,v,w*h))*h';
    sh=sum(h,2)';
%%% experiment: 
%     beta_w = 10000;
%     sh=sum(h,2)'+beta_w*sum(w.^2,1);
%%%
    w(:,peval.dovec_w)=bsxfun(@rdivide,w(:,peval.dovec_w).*(y2(:,peval.dovec_w)),sh(peval.dovec_w));
    w=max(w,eps); % adjust small values to avoid undeflow
        
    % normalization:
    sw = sum(w,1)'; % summation after the update of w...
    w = bsxfun(@rdivide, w, sw');
    h = bsxfun(@times, h, sw); % this is to keep the product WH constant after normalisation of w        
    if rem(ii,checkddivfreq)==0
        d = ddivergence(v, w*h);
        fprintf('Cycle %g D-divergence %g\n',ii-1,d)
        if peval.showprogress
            dall(indexd)=d;
            plotprogress(w,dall(1:indexd),peval)
            indexd=indexd+1;
            
        end
        if d < peval.ddterm
            break
        end
    end    
end
peval.ddiv_end = ddivergence(v, w*h); % final values of the d-divergence
peval.numiter = ii;
peval.maxiter_reached_flag = 0;
if ii == peval.maxiter
    mfprintf(peval.fid,'\nMAximum number of iteration (%g) reached! \n', peval.maxiter)
    peval.maxiter_reached_flag = 1;
end

varargout = struct('w',w,'h',h);
end % Main Function

% Nested functions: 
function peval=setDefaultValuesPeval(peval)

if ~isfield(peval, 'dterm'); peval.dterm = 1; end %termination criterion
if ~isfield(peval, 'ddterm'); peval.ddterm = 1; end %termination criterion
if ~isfield(peval, 'maxiter'); peval.maxiter = 1000; end
if ~isfield(peval,'showimage'); peval.showimage =0; end %showing progres images
if ~isfield(peval,'bgcomp'); peval.bgcomp = 1; end %last componnent as a backgroufn
if peval.bgcomp
    if ~isfield(peval, 'fix_bg_w') 
        peval.fix_bg_w=1; % last component w is fixed and not updated       
    end
    if ~isfield(peval, 'fix_bg_h')    
        peval.fix_bg_h=1; % last component h is fixed and not updated       
    end
end
end

function checkv(v)
% test for negative values in v
if min(v(:)) < 0
    error('matrix entries can not be negative');
end
if min(sum(v,2)) == 0
    error('not all entries in a row can be zero');
end
end

function [dovec_w, dovec_h]=setDoVec(peval)
dovec_w=1:peval.ncomp;
dovec_h=1:peval.ncomp;
if peval.fix_bg_w
    dovec_w=1:peval.ncomp-1;
end
if peval.fix_bg_h
    dovec_h=1:peval.ncomp-1;
end
if isfield(peval, 'fix_w')
    dovec_w = removerows(dovec_w',peval.fix_w)';
end
if isfield(peval, 'fix_h')
    dovec_h = removerows(dovec_h',peval.fix_h)';
end
if peval.fix_bg_w
    mfprintf(peval.fid, '''w'' fixed for background component [%g].\n',peval.ncomp)
end
if peval.fix_bg_h
     mfprintf(peval.fid, '''h'' fixed for background component [%g].\n',peval.ncomp)
end

mfprintf(peval.fid, '''w'' will be updated for components: ')
mfprintf(peval.fid, '['); mfprintf (peval.fid, '%g ', dovec_w); mfprintf (peval.fid, ']\n');                        
mfprintf(peval.fid, '''h'' will be updated for components: ')
mfprintf(peval.fid, '['); mfprintf (peval.fid, '%g ', dovec_h); mfprintf (peval.fid, ']\n');                        
end

function plotprogress(w,d,peval)
    imstiled(reshape(w, peval.nx, peval.ny, peval.ncomp),10, 'gray',[],[],1)
    plot(d(2:end))
end