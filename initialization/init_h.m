function h=init_h(method,peval,image,dpixc)
% h=init_h(method,peval,image,dpixc)
% method:   'rand' %random initialization
%           'image' %specified initialization for example image=double(array2im(dpixc_ind))       
%           'image_repmat' %for example image = mean(image,3);        
%           'res' %image is res.h -> initialisation from the nmf results
% addbackgroundcomponent:   1 adds one component flat component as a background (peval.ncomp th)
%                           0 no background component
if ~isfield(peval, 'bgcomp')
    peval.bgcomp = 1;
    mfprintf(peval.fid, 'Value of the peval.bgcomp was set to %g (by default the last component is backgroud).\n',peval.bgcomp)
end

switch method
    case 'rand'
        h = mean(dpixc(:)-peval.bg)*rand(peval.ncomp,peval.nt);       
        msg='h initialzied as uniform random.';
    case 'res'
        h = image;
        msg='h initialized from the results res.h.';        
end
h=max(h, eps); % To avoid zeros...

if isfield (peval,'fid')
    mfprintf(peval.fid, [msg '\n'])
else
    fprintf([msg '\n']);
end

if peval.bgcomp
    h(peval.ncomp,:)=peval.nx*peval.ny*peval.bg;
    mfprintf(peval.fid, 'Last component [%g] of h initialised as a flat background. (background=%g)\n',peval.ncomp,peval.bg);    
end
