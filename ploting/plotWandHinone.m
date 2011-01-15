ncomp = 5;
lsub = 10; %length of subtitles
for ii=1:ncomp
    indexim = (ii-1)*lsub+1; 
    subplot(ncomp,lsub, indexim)
    ims(double(dpixc_ind{ii}), 'gray', 0, 1);
    subplot(ncomp,lsub,[indexim+1, indexim+lsub-1]);
    plot(res.h(ii,1:200))
    set(gca,'xtick',[],'ytick',[])
end