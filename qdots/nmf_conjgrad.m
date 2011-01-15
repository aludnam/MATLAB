function varargout=nmf_conjgrad(v,winit,hinit,peval,verbose)
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
% w    : N x r NMF factor
% h    : r x M NMF factor
%
% winit - initial value for w
% hinit - initial value for h
% peval.h_fixvec & peval.w_fixvec:
% fixvec - vector which component should be fixed: eg [2 3] will
% fix second and third component while varying the first...
fprintf('Conjgrad NMF iterations\n')

if ~isfield(peval, 'maxh'); peval.maxh=10; end
if ~isfield(peval, 'maxw'); peval.maxw=10; end
if ~isfield(peval, 'ddterm'); peval.ddterm = 1; end
if ~isfield(peval, 'maxiter'); peval.maxiter = 100; end

options_h = zeros(1,18);
options_h(1)=1;                  %to display error values
options_h(7)=1;
options_h(9)=0;                   %to check gradient
options_h(14)=peval.maxh;         %maximum number of iterations

options_w = options_h;
options_w(14)=peval.maxw;

% test for negative values in v
if min(min(v)) < 0
    error('matrix entries can not be negative');
end
if min(sum(v,2)) == 0
    error('not all entries in a row can be zero');
end

ix0 = find(v==0);
if ~isempty(ix0)
    v(ix0) = eps; %othervise problems in Ddiv....
    fprintf('%g pixel values changed from 0 to %g\n', length(ix0), eps)
end

if ~isempty(peval.w_fixvec)
    fprintf('Fixing [ ');
    fprintf('%g ', peval.w_fixvec);
    fprintf('] component of ''w''\n');
end

if ~isempty(peval.h_fixvec)
    fprintf('Fixing [ ');
    fprintf('%g ', peval.h_fixvec);
    fprintf('] component of ''h''\n');
end

w_vec = 1:size(winit,2);
h_vec = 1:size(hinit,1);

peval.w_dovec = find(~(w_vec==peval.w_fixvec));
peval.h_dovec = find(~(h_vec==peval.h_fixvec));

[winitn, hinitn] = nomalize_wh(winit,hinit);
w = winitn(:, peval.w_dovec);
h = hinitn(peval.h_dovec, :);
w_fix = winitn(:, peval.w_fixvec);
h_fix = hinitn(peval.h_fixvec, :);

[n,m]=size(v);
kw = length(peval.w_dovec);
kh = length(peval.h_dovec);

d(1)=ddivergence(v, [w w_fix]*[h;h_fix]);
% d(1) = [];
dd(1)=realmax;
ii=2;
while and(dd>peval.ddterm, ii<peval.maxiter)
    if ~isempty(h)
        [rHkte, options_h, flog_h, pointlog_h] = conjgrad('ddivHexp', reshape(log(h),1,m*kh), options_h, 'gradddivHexp', v, w, w_fix, h_fix, peval);
        h = exp(reshape(rHkte,kh,m));
    end
    if ~isempty(w)        
        [rWxke, options_w, flog_w, pointlog_w] = conjgrad('ddivWexp', reshape(log(w),1,n*kw), options_w, 'gradddivWexp', v, h, w_fix, h_fix, peval);
        w = exp(reshape(rWxke,n,kw));
        [w, h] = nomalize_wh(w,h);
    end
    
    w_all = [w, w_fix];
    h_all = [h; h_fix];   
%     d(ii) = ddivergence(v, w_all * h_all);
% mess:
    d = [d flog_h'];
    jj = length(flog_h);
    dd(ii) = abs(d(jj)-d(jj-1));
    htrace=reshape(pointlog_h,numel(pointlog_h)/(m*kh),kh,m);
%     dd(ii) = abs(d(ii)-d(ii-1));
    if verbose
        fprintf('Cycle %g D-divergence %g\n',ii-1,d(ii))
    end
    if rem(ii,2)==0
%         pause(0.5)        
        wr = reshape(w,peval.nx, peval.ny, peval.ncomp);
        imstiled (wr,1,[],1)
% % %         fprintf('reseting H!\n')
% % %         sh = size(h);
% % %         mh = mean(h(:));
% % %         h = mh*2*rand(sh);
    end
    
    ii=ii+1;
end

varargout{1}=w_all;
varargout{2}=h_all;
varargout{3}=peval;
varargout{4}=d;
varargout{5}=htrace;
