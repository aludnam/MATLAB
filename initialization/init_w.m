function [w,winit_pix]=init_w(method,peval,image)
% [w,winit_pix]=init_w(method,peval,image)
% method:   'rand' %random initialization
%           'image' %specified initialization for example image=double(array2im(dpixc_ind))       
%           'image_repmat' %for example image = mean(image,3);        
%           'res' %image is res.w -> initialisation from the nmf results
% addbackgroundcomponent:   1 adds one component flat component as a background (peval.ncomp th)
%                           0 no background component
if ~isfield(peval, 'bgcomp')
    peval.bgcomp = 1;
end

switch method
    case 'rand'
        winit_pix = normalize(rand(peval.nx,peval.ny,peval.ncomp));        
        msg='W initialzied as uniform random.';
    case 'image'
        winit_pix = normalize(image);
        msg='W initialized with specified 3D image.';
    case 'image_repmat'
        winit_pix(:,:,1:peval.ncomp) = normalize(repmat(image,[1,1,peval.ncomp]));
        msg='W initialized as repmat(image,[1,1,peval.ncomp].)';
    case 'res'
        s2=size(image,2); %it is only 2D (#pixX#comp)
        winit_pix(:,:,1:s2) = reshape(image, peval.nx, peval.ny, s2);
        msg='W initialized from the results res.w.';
end
winit_pix=max(winit_pix, eps); % To avoid zeros...
sw=size(winit_pix);
w=reshape(winit_pix,sw(1)*sw(2), size(winit_pix,3)); % it must be specifically set to size(winit_pix,3) for the case tehre is only one compoenent (ncomp=1)

if isfield (peval,'fid')
    mfprintf(peval.fid, [msg '\n'])
else
    fprintf([msg '\n']);
end

if peval.bgcomp
    w(:,peval.ncomp)=normalize(ones(sw(1)*sw(2),1));
    mfprintf(peval.fid, 'Last component [%g] of w initialised as a flat background. (background=%g)\n',peval.ncomp,peval.bg);    
end