function [a,b] = init_ab(dpixc, alph, beta, peval)
% [a,b] = init_ab(dpixc, alph, beta, peval)

 
 b = repmat(beta'+1, 1, peval.nt);
 a = repmat((sum(alph) + sum(dpixc(:)))/(peval.ncomp*peval.nt), peval.ncomp, peval.nt);
 if peval.addbgcomp
        ncomp = peval.ncomp-1;
        a = repmat((sum(alph) + sum(dpixc(:)-peval.bg))/(ncomp*peval.nt), ncomp, peval.nt);
        a(peval.ncomp,:)=peval.bg*ones(1,peval.nt)*peval.nx*peval.ny;        
 end
    