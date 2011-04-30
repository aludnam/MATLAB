function [varargout, peval]=updates_variational(peval, varargin)
% [varargout, peval]=variationalupdates(peval, varargin)
% Computes variational updates
% dvec = varargin{1};   % data (N x T)
% w = varargin{2};      % initialization for basis (N x peval.ncomp)
% h = varargin{3};      % initialization for a (peval.ncomp x T) and b
% a = h.a;      % initialization for a 
% b = h.b;
mfprintf(peval.fid, 'Variational updates.\n')
peval=setDefaultValuesPeval(peval);
dvec = varargin{1};
w = varargin{2};        % initialization for basis (N x peval.ncomp)
h = varargin{3};        % initialization for a (peval.ncomp x T) and b
a = h.a;                % initialization for a 
b = h.b;

alph=repmat(peval.alpha,1,peval.ncomp); % prior parametes for Gamma distribution (peval.ncomp x T)
beta=repmat(peval.beta,1,peval.ncomp);  % prior parametes for Gamma distribution (peval.ncomp x T)

[dovec_w, dovec_a]=setDoVec(peval); %sets which components used in updates
lb=zeros(1,peval.maxiter);
lb(1)=sum(lowerbound(dvec, w, alph, beta, a, b));
for updateindex=1:peval.maxiter;
%     n=update_n(w, a, b);
    a(dovec_a,:)=update_anew(dvec,n(:,dovec_a,:),alph(:,dovec_a), a(dovec_a,:), b);
    w(:,dovec_w)=update_w(dvec, n(:,dovec_w,:));
    lb(updateindex+1)=sum(lowerbound(dvec, w, alph, beta, a, b));
    if peval.showprogress
        plotprogress(w,lb(2:updateindex),peval)
    end
end

varargout = struct('w',w,'a',a,'lb',lb(end));

end % Main function

% Nested fucntions:
function peval=setDefaultValuesPeval(peval)
if ~isfield(peval,'showprogress')
    peval.showprogress = 0;
end

if ~isfield(peval, 'fix_bg_w')
    peval.fix_bg_w=0;
    if peval.addbgcomp
        peval.fix_bg_w=1;
    end
end
if ~isfield(peval, 'fix_bg_a')
    peval.fix_bg_a=0;
    if peval.addbgcomp
        peval.fix_bg_a=1;
    end
end
if peval.fix_bg_w
    mfprintf(peval.fid, '''w'' fixed for background component [%g].\n',peval.ncomp)
end
if peval.fix_bg_a
    mfprintf(peval.fid, '''a'' fixed for background component [%g].\n',peval.ncomp)
end
end

function [dovec_w, dovec_a]=setDoVec(peval)
dovec_w=1:peval.ncomp;
dovec_a=1:peval.ncomp;
if peval.fix_bg_w
    dovec_w=1:peval.ncomp-1;
end
if peval.fix_bg_a
    dovec_a=1:peval.ncomp-1;
end
if isfield(peval, 'fix_w')
    dovec_w = removerows(dovec_w',peval.fix_w)';
end
if isfield(peval, 'fix_a')
    dovec_a = removerows(dovec_a',peval.fix_a)';
end
mfprintf(peval.fid, '''w'' will be updated for components: ')
mfprintf(peval.fid, '['); mfprintf (peval.fid, '%g ', dovec_w); mfprintf (peval.fid, ']\n');
mfprintf(peval.fid, '''a'' will be updated for components: ')
mfprintf(peval.fid, '['); mfprintf (peval.fid, '%g ', dovec_a); mfprintf (peval.fid, ']\n');
end

function plotprogress(w,lb,peval)
imstiled(reshape(w, peval.nx, peval.ny, peval.ncomp),10, 'gray',[],[],1)
plot(lb)
end