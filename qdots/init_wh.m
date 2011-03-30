function [w,h,winit_pix]=init_wh(methodw,methodh,peval,imagew,imageh,dpixc,betain,alphain)

updates = func2str(peval.fun);

switch updates
    case {'updates_nmfclassic', 'updates_variational'}        
        [w,winit_pix]=init_w(methodw,peval,imagew);
        h=init_h(methodh,peval,imageh);
    case 'updates_map'
        w=init_wmap(method,peval,imagew);
        winit_pix=makegauss(reshape(w,1,2*peval.ncomp),peval.sigmaPSF, [peval.nx, peval.ny]);
        h = init_ab(dpixc, alphain, betain, peval);
        
end