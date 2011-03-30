function [w,winit_pix]=init_w_general(method,peval,image)

updates = func2str(peval.fun);

switch updates
    case {'updates_nmfclassic', 'updates_variational'}        
        [w,winit_pix]=init_w(method,peval,image);
    case 'updates_map'
        w=init_wmap(method,peval,image);
        winit_pix=makegauss(reshape(w,1,2*peval.ncomp),peval.sigmaPSF, [peval.nx, peval.ny]);
end