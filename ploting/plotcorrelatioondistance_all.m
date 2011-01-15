function maxcorrel=plotcorrelatioondistance_all(delta_vec,dir1,dir2, name1, name2, marker, col, savethis, correlationmethod, niter)
% maxcorrel=plotcorrelatioondistance(delta_vec,res_dir, res_name, marker,
% color savethis)
if ~exist('correlationmethod', 'var')
    correlationmethod = 'pos'; % positive correlations
end
if ~exist('niter','var')
    niter=10;
end

method='average';
ii=1;
for delta = delta_vec
    
    for iteration =1:niter;
        res_dir = [dir1 num2str(delta*100) dir2];
        res_name = [name1 num2str(iteration) name2];
        load ([res_dir '/' res_name]);
%         load ([peval.data_path peval.data_dir '/' peval.data_name]); %data
        [Wxk,Hkt,centers,Vxkpix]=reshapeGaP(res.hvec,res.cxcy,peval);
        Vxtpixbg=reshape(Wxk*Hkt,peval.nx,peval.ny,peval.nt)+peval.bg;
        resid=(Vxtpixbg-res.dpixc);
        resid_norm=resid./sqrt(Vxtpixbg);
        %         [Z,H,T,perm] = dendrogram_subtreepixels(resid_norm,'average', p , centers(:,1)+1, centers(:,2)+1, savethis)
        data=resid_norm;
        sized = size(data);
        dveccr= reshape(data,sized(1)*sized(2), sized(3));
        ccd = (corrcoef(dveccr'));
        ccds = correlation2distance(ccd, correlationmethod);
%         ccds=squareform(1-abs(ccd));
%         e=eye(size(ccd));
%         ccdflip=e-(ccd-e);
%         ccds=squareform(1-ccdflip);
%         
        Z = linkage(ccds,method);
        
        if savethis
            save ([res_dir '/corrcoefdist.mat'], 'Z')
        end
        maxcorrel(ii,iteration)=1-min(Z(:,3));
    end
    ii=ii+1;
end


% style='o--r';
% marker='o';
% col='r';
style = [marker col '--'];
figure(1)
hold on
deltanm=delta_vec*106;
m=mean(maxcorrel,2);
s=std(maxcorrel,[],2);
mincorr=min(maxcorrel,[],2);
errorbar(deltanm, m,s,style)
plot(deltanm,mincorr,['-' marker col], 'linewidth',2)
xlabel('Separation of sources [nm]')
ylabel('Maximum correlation in residuals')
grid on


figure(2)
hold on
m=mean(maxcorrel,2);
s=std(maxcorrel,[],2);
errorbar(delta_vec, m,s,style)
xlabel('Separation of sources [pix]')
plot(delta_vec,mincorr,['-' marker col], 'linewidth',2)
ylabel('Maximum correlation in residuals')
grid on

    