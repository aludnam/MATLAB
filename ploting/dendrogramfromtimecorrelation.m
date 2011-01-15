function [Z,H,T,perm] = dendrogramfromtimecorrelation(data, method, savethis, name_figure)
% [Z,H,T,perm] = dendrogramfromtimecorrelation(data, method, savethis,name_figure)


if ~exist('savethis', 'var')
    savethis = 0;
    if ~exist('name_figure', 'var')
        name_figure = ['dendrogram_' method];
    end
    
end


dveccr= reshape(data,size(data,1)*size(data,2),size(data,3));
ccd = (corrcoef(dveccr'));
ccds=squareform(1-ccd);
Z = linkage(ccds,method);
[H,T,perm] = dendrogram(Z,0);
ylim([0 max(Z(:,3))])
set(H,'color','k')
grid on
if savethis
    SaveImageFULL(name_figure)
end
    
