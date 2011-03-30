function [varargout, peval]=updates_nmfclassic(peval,varargin)
% [varargout, peval]=updates_nmfclassic(peval,varargin)
%
% Jean-Philippe Brunet
% Cancer Genomics
% The Broad Institute
% brunet@broad.mit.edu
% NMF divergence update equations :
% Lee, D..D., and Seung, H.S., (2001), 'Algorithms for Non-negative Matrix
% Factorization', Adv. Neural Info. Proc. Syst. 13, 556-562.
%
% v (n,m) : N (genes) x M (samples) original matrix
%           Numerical data only.
%           Must be non negative.
%           Not all entries in a row can be 0. If so, add a small constant to the
%           matrix, eg.v+0.01*min(min(v)),and restart.
%
% r       : number of desired factors (rank of the factorization)
%
% verbose : prints iteration count and changes in connectivity matrix elements
%           unless verbose is 0
%
% Note : NMF iterations stop when connectivity matrix has not changed
%        for 10*stopconv interations. This is experimental and can be
%        adjusted.
%
% w    : N x r NMF factor
% h    : r x M NMF factor
%
% v = varargin{1};
% w = varargin{2};
% h = varargin{3};
% peval.h_fixvec & peval.w_fixvec:
% fixvec - vector which component should be fixed: eg [2 3] will
% fix second and third component while varying the first...

mfprintf(peval.fid,'Classic NMF iterations\n')
peval=setDefaultValuesPeval(peval);
v = varargin{1};
w = varargin{2};
h = varargin{3};
[peval.dovec_w, peval.dovec_h]=setDoVec(peval);

checkv(v) %check values of v
[n,m]=size(v);
checkddivfreq=10; % How often ot check ddiv
dall=zeros(1, ceil(peval.maxiter/checkddivfreq));
indexd=1;
for ii=2:peval.maxiter
    w_old = w;
    h_old = h;
        
    x1=repmat(sum(w,1)',1,m);
    y1=w'*(v./(w*h));
    h(peval.dovec_h,:)=h(peval.dovec_h,:).*(y1(peval.dovec_h,:))./x1(peval.dovec_h,:);
    h=max(h,eps); % adjust small values to avoid undeflow
    
    x2=repmat(sum(h,2)',n,1);
    y2=(v./(w*h))*h';
    w(:,peval.dovec_w)=w(:,peval.dovec_w).*(y2(:,peval.dovec_w))./x2(:,peval.dovec_w);
    w=max(w,eps); % adjust small values to avoid undeflow
        
    % normalization of all h:
    sumw = sum(w,1);
    w = w./repmat(sumw,n,1); %normalization of each component
    h = h.*repmat(sumw',1,m); %to keep the multiplication equal
        
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
if ~isfield(peval,'addbgcomp'); peval.addbgcomp = 1; end %last componnent as a backgroufn
if peval.bgcomp
    if ~isfield(peval, 'fix_bg_w')    
        peval.fix_bg_w=1;        
    end
    if ~isfield(peval, 'fix_bg_h')    
        peval.fix_bg_h=1;        
    end
end
end

function checkv(v)
% test for negative values in v
if min(min(v)) < 0
    error('matrix entries can not be negative');
    return
end
if min(sum(v,2)) == 0
    error('not all entries in a row can be zero');
    return
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
     mfprintf(peval.fid, '''g'' fixed for background component [%g].\n',peval.ncomp)
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