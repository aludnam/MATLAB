function varargout=nmf_canny(v,winit,hinit,peval,verbose)
% NMF updates from [Canny et all 2004]
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
% winit - initial value for w
% hinit - initial value for h
% peval.h_fixvec & peval.w_fixvec:
% fixvec - vector which component should be fixed: eg [2 3] will
% fix second and third component while varying the first...
fprintf('Classic NMF iterations\n')

if ~isfield(peval, 'dterm'); peval.dterm = 1; end %termination criterion
if ~isfield(peval, 'ddterm'); peval.ddterm = 1; end %termination criterion
if ~isfield(peval, 'maxiter'); peval.maxiter = 1000; end


% test for negative values in v
if min(min(v)) < 0
    error('matrix entries can not be negative');
    return
end
if min(sum(v,2)) == 0
    error('not all entries in a row can be zero');
    return
end

[n,m]=size(v);

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

w = winit;
h = hinit;

d(1) = ddivergence(v, w*h);
htrace(1,:,:) = hinit;

w_vec = 1:size(w,2);
h_vec = 1:size(h,1);
% peval.w_dovec = find(~(w_vec==peval.w_fixvec));
% peval.h_dovec = find(~(h_vec==peval.h_fixvec));
peval.w_dovec = find(~ismember(w_vec,peval.w_fixvec));
peval.h_dovec = find(~ismember(h_vec,peval.h_fixvec));
% f1=figure('outerposition',[10 445 560 420],'name','W');
% f2=figure('outerposition',[10 10 560 420],'name','objective function');
% f3=figure('outerposition',[570 10 560 420],'name','derivative objective function');

for ii=2:peval.maxiter
    
    w_old = w;
    h_old = h;
        
    x1=repmat(sum(w,1)',1,m)+1/peval.betaH;
    y1=w'*(v./(w*h)) + (peval.alphaH-1)./h;
    h(peval.h_dovec,:)=h(peval.h_dovec,:).*(y1(peval.h_dovec,:))./x1(peval.h_dovec,:);
    h=max(h,eps); % adjust small values to avoid undeflow
    
    x2=repmat(sum(h,2)',n,1);
    y2=(v./(w*h))*h';
    w(:,peval.w_dovec)=w(:,peval.w_dovec).*(y2(:,peval.w_dovec))./x2(:,peval.w_dovec);
    w=max(w,eps); % adjust small values to avoid undeflow
        
    % normalization of all h:
    sumw = sum(w,1);
    w = w./repmat(sumw,n,1); %normalization of each component
    h = h.*repmat(sumw',1,m); %to keep the multiplication equal
    
    if verbose
    if rem(ii,10)==0
%         pause(0.5)        
        wr = reshape(w,peval.nx, peval.ny, peval.ncomp+1);
%         figure(f1)
        imstiled (wr,1,[],1)        
%         [AX,H1,H2] =plotyy([1:ii-1], log(d./max(d)), [1:ii-1], log(dd./max(dd)));
        [AX,H1,H2] =plotyy([1:ii-1], log(d./max(d)), [1:ii-1], log(wdiff./max(wdiff)));
        set(AX(1),'xtick',[],'ytick',[])
        set(AX(2),'xtick',[],'ytick',[])
        set(H2,'LineStyle',':')
        set(get(AX(1),'Ylabel'),'String','d')
        set(get(AX(2),'Ylabel'),'String','w diff')

%         xlim([min(50,ii-100), ii])
        
%         figure(f2)
%         semilogy(d)
%         figure(f3)
%         semilogy(dd)
% % %         fprintf('reseting H!\n')
% % %         sh = size(h);
% % %         mh = mean(h(:));
% % %         h = mh*2*rand(sh);
        fprintf('Cycle %g D-divergence %g\n',ii-1,d(ii-1))
    end
    end
    d(ii) = ddivergence(v, w*h);
    dd(ii) = abs(d(ii)-d(ii-1));
    wdiff(ii) = sum((w_old(:)-w(:)).^2);
    htrace(ii,:,:)=h;
%     wtrace(:,:,ii) = w;
%     htrace(:,:,ii) = h;
     
    if dd(ii) < peval.ddterm
        break
    end
    if d(ii) < peval.dterm
        break
    end    
%     if verbose
%         fprintf('Cycle %g D-divergence %g\n',ii-1,d(ii))
%     end
end

peval.numiter = ii;
peval.maxiter_reached_flag = 0;
if ii == peval.maxiter
    fprintf('\nMAximum number of iteration (%g) reached! \n', peval.maxiter)
    peval.maxiter_reached_flag = 1;
end

varargout{1}=w;
varargout{2}=h;
varargout{3}=peval;
varargout{4}=d;
varargout{5}=htrace;