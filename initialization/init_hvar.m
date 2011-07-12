function h = init_hvar(method, peval, image, dpixc)
% h = init_hvar(method, peval, image, dpixc)
% Initaialization of the a, b parameters of the factorised distributions in
% the variational updates of hte GaP model.
% method: 'res' initialization from the specific values in image (eg
% image=res.h)
% dpixc: data

sumdata = sum(dpixc(:));
alph=repmat(peval.alpha,1,peval.ncomp); % prior parametes for Gamma distribution (peval.ncomp x T)
beta=repmat(peval.beta,1,peval.ncomp);  % prior parametes for Gamma distribution (peval.ncomp x T)

b = repmat(beta'+1, 1, peval.nt);
switch method
    case 'res'
        a = image;
        msg='a initialized from the results res.h.';
    case 'image'
        a = image;
        msg='a initialized from the specified image.';
    otherwise
        a = repmat((sum(alph) + sumdata)/(peval.ncomp*peval.nt), peval.ncomp, peval.nt);
        msg='a initialized uniformly constant';
        if peval.bgcomp
            ncomp = peval.ncomp-1;
            a = repmat((sum(alph) + sum(sumdata-peval.bg))/(ncomp*peval.nt), ncomp, peval.nt);
            a(peval.ncomp,:)=peval.bg*ones(1,peval.nt)*peval.nx*peval.ny;
            mfprintf(peval.fid, 'Last component [%g] of h initialised as a flat background. (background=%g)\n',peval.ncomp,peval.bg);
        end
end
if isfield (peval,'fid')
    mfprintf(peval.fid, [msg '\n'])
else
    fprintf([msg '\n']);
end



h = struct('a',a,'b', b);