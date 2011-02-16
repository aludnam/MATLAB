function [a,b] = init_ab(dpixc, alph, beta, peval)
% [a,b] = init_ab(dpixc, alph, beta, peval)

 a = repmat((sum(alph) + sum(dpixc(:)))/(peval.ncomp*peval.nt), peval.ncomp-1, peval.nt);
 b = repmat(beta'+1, 1, peval.nt);
 if peval.addbgcomp
        a(peval.ncomp,:)=peval.bg*ones(1,peval.nt)*peval.nx*peval.ny;        
 end
    