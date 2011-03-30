function w=init_wmap(method,peval,image)
% w=init_wmap(method,peval,image)
% method:   'rand' %random initialization
switch method
    case 'rand'
        w = [peval.nx*rand(peval.ncomp, 1), peval.ny*rand(peval.ncomp,1)];
        msg='W initialzied as uniform random in the image.';
end


fprintf([msg '\n']);


