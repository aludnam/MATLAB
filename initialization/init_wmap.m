function w=init_wmap(method,peval,image)
% w=init_wmap(method,peval,image)
% method:   'rand' %random initialization
%           'res' initialization fro the values of image
switch method
    case 'rand'
        w = [peval.nx*rand(peval.ncomp, 1), peval.ny*rand(peval.ncomp,1)];
        msg='W initialzied as uniform random in the image.';
    case 'res'
        w = image;
        msg='W initialized from the supplied values.';        
end

mfprintf(peval.fid,[msg '\n']);


