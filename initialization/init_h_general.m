function h=init_h_general(method,peval,image,dpixc)
% h=init_h_general(method,peval,image,dpixc)
% Initialization of the h (NMF, variationla) or a (h.a) and b (h.b) for MAP
% updates. 
% NMF and MAP:
% method:   'rand' %random initialization
%           'image' %specified initialization for example image=double(array2im(dpixc_ind))       
%           'image_repmat' %for example image = mean(image,3);        
%           'res' %image is res.w -> initialisation from the nmf results
% addbackgroundcomponent:   1 adds one component flat component as a background (peval.ncomp th)
%                           0 no background component
% Variational: 
% method:   not used
% image:    sum(dpixc(:)) (sum of the data)

updates = func2str(peval.fun);

switch lower(updates)
    case {'updates_nmfclassic', 'updates_map'}        
        h=init_h(method,peval,image, dpixc);
    case 'updates_variational'
        h=init_hvar(method, peval, image, dpixc); %here image is sum(dpixc(:))
end