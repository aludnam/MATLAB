function ccds = correlation2distance(ccd, method)
%ccds = correlation2distance(ccd)
%   Computes distance matrix from the correlation coefficients matrix (ccd = (corrcoef(data));)
%   method :    'pos' %only positive correlations
%               'neg' %negative correlations (anitcorrelations)
%               'abs' %abs value of correlation coefficients
%               'sqr' %square value of correlation coefficients

switch method
    case 'pos' %positive correlations
        ccds=squareform(1-ccd);
        fprintf('Positive correlations coefficients.\n')
    case 'neg' %negative correlations (anitcorrelations)
        e=eye(size(ccd));
        ccdflip=e-(ccd-e);
        ccds=squareform(1-ccdflip);
        fprintf('Negative correlations coefficients.\n')
    case 'abs' %abs value of correlation coefficients
        ccds=squareform(1-abs(ccd));
        fprintf('Absolute values of correlations coefficients.\n')
    case 'sqr' %square of correlation coefficients       
        ccds=squareform(1-ccd.^2);
        fprintf('Square values of correlations coefficients.\n')
end