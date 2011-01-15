function [w,h, wtrace, htrace, d]=nmf_test(v,r,verbose,winit,hinit,w_fixvec,h_fixvec)
%
% Jean-Philippe Brunet
% Cancer Genomics
% The Broad Institute
% brunet@broad.mit.edu
%
% This software and its documentation are copyright 2004 by the
% Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
% This software is supplied without any warranty or guaranteed support whatsoever.
% Neither the Broad Institute nor MIT can not be responsible for its use, misuse,
% or functionality.
%
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
% winit - initial value for w
% hinit - initial value for h
% h_fixvec & w_fixvec:
% fixvec - binary vector which component should be fixed: eg [2 3] will
% fix second and third component while varying the first...

if ~exist ('w_fixvec', 'var')
    w_fixvec = [];
end

if ~exist ('h_fixvec', 'var')
    h_fixvec = [];
end


% test for negative values in v
if min(min(v)) < 0
    error('matrix entries can not be negative');
    return
end
if min(sum(v,2)) == 0
    error('not all entries in a row can be zero');
    return
end

ix0 = find(v==0);
v(ix0) = eps; %othervise problems in Ddiv....
     

[n,m]=size(v);
stopconv=40;      % stopping criterion (can be adjusted)
niter = 1000;     % maximum number of iterations (can be adjusted)

cons=zeros(m,m);
consold=cons;
inc=0;
j=0;

%
% initialize random w and h
%
w = rand(n,r);
if exist('winit', 'var')
    if ~isempty(winit)
        w = winit;
        fprintf ('initial values of ''w'' used...\n');
    end
end

h = rand(r,m);
if exist('hinit', 'var')
    if ~isempty(hinit)
        hinit = hinit./repmat(sum(hinit,2),1,m); %normalization of initial components
        h = hinit;
        fprintf ('initial values of ''h'' used...\n');
    end    
end

if exist('w_fixvec', 'var')
    fprintf('Fixing [ ');
    fprintf('%g ', w_fixvec);
    fprintf('] component of ''w''\n');
end

if exist('h_fixvec', 'var')
    fprintf('Fixing [ ');
    fprintf('%g ', h_fixvec);
    fprintf('] component of ''h''\n');
end

% oldpar = load('/afs/inf.ed.ac.uk/user/s08/s0880377/project/MATLAB/qdots/h5.mat');
% w = oldpar.w;
% h=oldpar.h;

h = h./repmat(sum(h,2),1,m); %normalization of each component

htrace(:,:,1) = h;
wtrace(:,:,1) = w;
d(1) = Ddiv(v, w*h);
dterm = 1; %can be changed

w_vec = 1:size(w,2);
w_dovec = find(~(w_vec==w_fixvec));

h_vec = 1:size(h,1);
h_dovec = find(~(h_vec==h_fixvec));


for i=2:niter
    
    w_old = w;
    h_old = h;
    
    % divergence-reducing NMF iterations
    
%     x2=repmat(sum(h,2)',n,1);
%     w=w.*((v./(w*h))*h')./x2;
%     w=max(w,eps); % adjust small values to avoid undeflow
%     w(:,w_fixvec) = winit(:,w_fixvec);
%     wtrace(:,:,i) = w;

%     x1=repmat(sum(w,1)',1,m);
%     h=h.*(w'*(v./(w*h)))./x1;
%     h=max(h,eps); % adjust small values to avoid undeflow
%     h(h_fixvec,:) = hinit(h_fixvec,:); 
%     htrace(:,:,i) = h;



    x2=repmat(sum(h,2)',n,1);
    y2=(v./(w*h))*h';
    w(:,w_dovec)=w(:,w_dovec).*(y2(:,w_dovec))./x2(:,w_dovec);
    w=max(w,eps); % adjust small values to avoid undeflow
    wtrace(:,:,i) = w;
    
% % %     x1=repmat(sum(w,1)',1,m);
% % %     y1=w'*(v./(w*h));
% % %     h(h_dovec,:)=h(h_dovec,:).*(y1(h_dovec,:))./x1(h_dovec,:);
% % %     h=max(h,eps); % adjust small values to avoid undeflow
% % %     
% % %     % normalization of all h:
% % %     sumh = sum(h,2);
% % %     h = h./repmat(sumh,1,m); %normalization of each component
% % %     w = w.*repmat(sumh',n,1); %to keep the multiplication equal
    
    htrace(:,:,i) = h;
    
    d(i) = ddivergence(v, w*h);
    dd(i) = abs(d(i)-d(i-1));
    if dd(i) < dterm
        break
    end
    
    
%     wterm=mean(winit(:))/(10^3);
%     hterm=mean(hinit(:))/(10^3);
%     dh(i) = max(abs(h_old(:)-h(:)));
%     dw(i) = max(abs(w_old(:)-w(:)));

    %     
    %     fprintf('maxDh %g\n',dh(i))
    %     fprintf('maxDw %g\n',dw(i))
    fprintf('Ddiv %g\n',d(i))
%     if dh(i) < hterm &  dw(i) < wterm
%         break
%     end
    

    % test convergence every 10 iterations
%     if(mod(i,10)==0)
%         j=j+1;
%         htrace(:,:,j+1) = h;
%         wtrace(:,:,j+1) = w;
%         
%         % construct connectivity matrix
%         [y,index]=max(h,[],1);   %find largest factor
%         mat1=repmat(index,m,1);  % spread index down
%         mat2=repmat(index',1,m); % spread index right
%         cons=mat1==mat2;
%         
%         if(sum(sum(cons~=consold))==0) % connectivity matrix has not changed
%             inc=inc+1;                     %accumulate count
%         else
%             inc=0;                         % else restart count
%         end
%         if verbose                     % prints number of changing elements
%             fprintf('\t%d\t%d\t%d\n',i,inc,sum(sum(cons~=consold))),
%         end
%         
%         if(inc>stopconv)
%             break,                % assume convergence is connectivity stops changing
%         end
%         
%         consold=cons;
%         
%     end
end

if i == niter
    fprintf('\nMAximum number of iteration (%g) reached! \n', niter)
end
end