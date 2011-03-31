function [w,winit_pix]=init_w_general(method,peval,image)
% [w,winit_pix]=init_w_general(method,peval,image)
% Initialization of w (NMF, Variatioanl) or cooridinates of the gaussians
% (MAP)
% NMF & Variational:
% method:   'rand' %random initialization
%           'image' %specified initialization for example image=double(array2im(dpixc_ind))       
%           'image_repmat' %for example image = mean(image,3);        
%           'res' %image is res.w -> initialisation from the nmf results
% addbackgroundcomponent:   1 adds one component flat component as a background (peval.ncomp th)
%                           0 no background component
% MAP:
% method:   'rand' ranodm cooridnates in the images


updates = func2str(peval.fun);

switch lower(updates)
    case {'updates_nmfclassic', 'updates_variational'}        
        [w,winit_pix]=init_w(method,peval,image);
    case 'updates_map'
        w=init_wmap(method,peval,image);
        winit_pix=makegauss(reshape(w,1,2*peval.ncomp),peval.sigmaPSF, [peval.nx, peval.ny]);
end