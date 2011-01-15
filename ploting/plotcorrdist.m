function plotcorrdist(maxcorrel, delta_vec,col, marker)
%plotcorrdist(maxcorrel, delta_vec,col, marker)

style = [marker col '--'];
figure(1)
hold on
deltanm=delta_vec*106;
m=mean(maxcorrel,2);
s=std(maxcorrel,[],2);
mincorr=min(maxcorrel,[],2);
errorbar(deltanm, m,s,style)
plot(deltanm,mincorr,['-s' col], 'linewidth',2)
xlabel('Separation of sources [nm]')
ylabel('Maximum correlation in residuals')
grid on


figure(2)
hold on
m=mean(maxcorrel,2);
s=std(maxcorrel,[],2);
errorbar(delta_vec, m,s,style)
xlabel('Separation of sources [pix]')
plot(delta_vec,mincorr,['-s' col], 'linewidth',2)
ylabel('Maximum correlation in residuals')
grid on