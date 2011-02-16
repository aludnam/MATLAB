function [varargout, peval]=variationalupdates(peval, varargin)
% [varargout, peval]=variationalupdates(peval, varargin)
% Computes variational updates
% dvec = varargin{1};   % data (N x T)
% w = varargin{2};      % initialization for basis (N x peval.ncomp)
% alph = varargin{3};   % prior parametes for Gamma distribution (peval.ncomp x T)
% beta = varargin{4};   % prior parametes for Gamma distribution
% a = varargin{5};      % initialization for a (peval.ncomp x T)
% b = varargin{6};
if ~isfield(peval,'showprogress')
    showprogress = 0;
else 
    % plots the progress of the evaluation
    showprogress = peval.showprogress;
end

dvec = varargin{1};
w = varargin{2};
alph = varargin{3};
beta = varargin{4};
a = varargin{5};
b = varargin{6};


lb(1)=sum(lowerbound(dvec, w, alph, beta, a, b));
for updateindex=1:peval.maxiter;
    n=update_n(w, a, b);
    a(1:peval.ncomp-1,:)=update_a(dvec,n(:,1:peval.ncomp-1,:),alph(:,1:peval.ncomp-1));
    w(:,1:peval.ncomp-1)=update_w(dvec, n(:,1:peval.ncomp-1,:));
    lb(updateindex+1)=sum(lowerbound(dvec, w, alph, beta, a, b));
    if showprogress
        imstiled(reshape(w, peval.nx, peval.ny, peval.ncomp),10, 'gray',[],[],1)
        plot(lb(2:end))
    end
end

varargout = struct('w',w,'a',a,'lb',lb);