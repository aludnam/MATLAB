function [w,winit_pix]=init_w(method,peval,image)
% [w,winit_pix]=init_w(method,peval,image)
% method:   'rand' %random initialization
%           'image' %specified initialization for example image=double(array2im(dpixc_ind))       
%           'image_repmat' %for example image = mean(image,3);        
%           'res' %image is res.w -> initialisation from the nmf results

switch method
    case 'rand'
        winit_pix = normalize(rand(peval.nx,peval.ny,peval.ncomp));        
        msg='W initialzied as uniform random.';
    case 'image'
        winit_pix = normalize(image);
        msg='W initialized with specified 3D image.';
    case 'image_repmat'
        an=normalize(image);
        winit_pix(:,:,1:peval.ncomp) = repmat(an,[1,1,peval.ncomp]);
        msg='W initialized as repmat(image,[1,1,peval.ncomp].)';
    case 'res'
        winit_pix(:,:,1:peval.ncomp) = reshape(image(:,1:peval.ncomp), peval.nx, peval.ny, peval.ncomp);
        msg='W initialized from the results res.w.';
end
winit_pix=max(winit_pix, eps); % To avoid zeros...
w=reshape(winit_pix,peval.nx*peval.ny, size(winit_pix,3));

if isfield (peval,'fid')
    mfprintf(peval.fid, [msg '\n'])
else
    fprintf([msg '\n']);
end
