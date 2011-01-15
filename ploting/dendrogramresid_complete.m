% plots dendrogram and after selecting nodes, plot them and corresponding
% pixels in different coulours
ca
savethis = 0;
% method = 'average'
method = 'complete'
[Z,H,T,perm] = dendrogramfromtimecorrelation(residual,method); %'average');
% Z = dendrogramfromtimecorrelation(residual,'average');

clear('cl')
cl(1).ixZ = find(round(10^4*Z(:,3))==8807);
% cl(2).ixZ = find(round(10^4*Z(:,3))==8939);
% cl(3).ixZ = find(round(10^4*Z(:,3))==3483);
% cl(4).ixZ = find(round(10^4*Z(:,3))==4619);
% cl(5).ixZ = find(round(10^4*Z(:,3))==6097);%1804);

% [H,T,perm] = dendrogram(Z,0);
co= get(gca, 'colororder');
set(H,'color','k')

imz2 = zeros([peval.nx*peval.ny, 3]);
for ii=1:length(cl)
    [cl(ii).ixZ_vec, cl(ii).endleaves_vec] = recursivesubtree(Z, cl(ii).ixZ, [], []);
    set(H(cl(ii).ixZ_vec), 'color', co(ii,:))
    imz2(cl(ii).endleaves_vec,1)=co(ii,1);
    imz2(cl(ii).endleaves_vec,2)=co(ii,2);
    imz2(cl(ii).endleaves_vec,3)=co(ii,3);
end

if savethis
SaveImageFULL(['dendrogramres_' method '_nc' num2str(peval.ncomp)])
end


joinchannels('rgb',reshape(imz2,peval.nx,peval.ny,3))
hold on 
scatter(p.x_vec-0.5, p.y_vec-0.5,200,'xr')
scatter(x_mu-1, y_mu-1,'b')
if savethis
SaveImageFULL(['res-lefttoright_ ' method '_nc' num2str(peval.ncomp)])
end