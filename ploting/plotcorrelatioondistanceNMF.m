function [maxcorrel, maxcorrelhmax, maxcorrelhmin]=plotcorrelatioondistanceNMF(delta_vec,dir1,dir2, name1, name2, marker, col, savethis)
% [maxcorrel, maxcorrelhmx, maxcorrelhmin]=plotcorrelatioondistance(delta_vec,res_dir, res_name, marker,
% color savethis)

method='average';
ii=1;
for delta = delta_vec
    
    for iteration =1:10;
        res_dir = [dir1 num2str(delta*100) dir2];
        res_name = [name1 num2str(iteration) name2];
        load ([res_dir '/' res_name]);        
        resid=(res.w*res.h-reshape(res.dpixc, peval.nx*peval.ny, peval.nt));
        dveccr=resid./sqrt(res.w*res.h);                
        ccd = (corrcoef(dveccr'));
        ccds=squareform(1-ccd);
        
        Z = linkage(ccds,method);
        
        if savethis
            save ([res_dir '/corrcoefdist.mat'], 'Z')
        end
        maxcorrel(ii,iteration)=max(max(ccd-eye(size(ccd))));
        maxcorrelhmax(ii,iteration)= max(max(corr(res.h(1:end-1,:)')-eye(size(res.h,1)-1)));
        maxcorrelhmin(ii,iteration)= min(min(corr(res.h(1:end-1,:)')));
    end
    ii=ii+1;
end

plotcorrdist(maxcorrel, delta_vec,col,marker)

figure(3)
hold on
deltanm=delta_vec*106;
m=mean(maxcorrelhmax,2);
s=std(maxcorrelhmax,[],2);
errorbar(deltanm, m,s,'x-','color', col)
% plot(deltanm,mincorr,['-s' col], 'linewidth',2)
xlabel('d [nm]')
ylabel('Maximum correlation intensities (res.h)')
grid on

figure(4)
hold on
deltanm=delta_vec*106;
m=mean(maxcorrelhmin,2);
s=std(maxcorrelhmin,[],2);
errorbar(deltanm, m,s,'x-','color', col)
% plot(deltanm,mincorr,['-s' col], 'linewidth',2)
xlabel('d [nm]')
ylabel('Minimum correlation intensities (res.h)')
grid on

% % style='o--r';
% % marker='o';
% % col='r';
% style = [marker col '--'];
% figure(1)
% hold on
% deltanm=delta_vec*106;
% m=mean(maxcorrel,2);
% s=std(maxcorrel,[],2);
% mincorr=min(maxcorrel,[],2);
% errorbar(deltanm, m,s,style)
% plot(deltanm,mincorr,['-s' col], 'linewidth',2)
% xlabel('Separation of sources [nm]')
% ylabel('Maximum correlation in residuals')
% grid on
% 
% 
% figure(2)
% hold on
% m=mean(maxcorrel,2);
% s=std(maxcorrel,[],2);
% errorbar(delta_vec, m,s,style)
% xlabel('Separation of sources [pix]')
% plot(delta_vec,mincorr,['-s' col], 'linewidth',2)
% ylabel('Maximum correlation in residuals')
% grid on
% 
