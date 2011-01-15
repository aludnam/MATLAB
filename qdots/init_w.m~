function winit_pix=init_w(method,peval,image)
% winit_pix=init_w(method,peval,image)

switch method
    case 'rand' %random initialization
        winit_pix = rand(peval.nx,peval.ny,peval.ncomp); 
        % need to be normalized!
    case 'image_repmat' %for example image = mean(image,3);        
        an=normalize(image); %2d image
        winit_pix(:,:,1:peval.ncomp) = repmat(an,[1,1,peval.ncomp]);
    case 'res' %image is res.w -> initialisation from the nmf results
        winit_pix(:,:,1:peval.ncomp) = reshape(image(:,1:peval.ncomp), peval.nx, peval.ny, peval.ncomp);
end