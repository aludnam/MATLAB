function sumsq = sumsqcomp(dpixc, res, peval)
ROIx=peval.ROIx(1):peval.ROIx(2);
ROIy=peval.ROIy(1):peval.ROIy(2);
ROIz=peval.ROIz(1):peval.ROIz(2);
dpixcs = dpixc(ROIx,ROIy,ROIz); 
dveccs = reshape(dpixcs, peval.numpix, peval.nt);
htr =res.htrace;
for ii=1:size(htr,1)
    htr(ii,3,:)=res.h(3,:);
    itmp=res.w*squeeze(htr(ii,:,:));
    isq=(dveccs-itmp).^2;
    sumsq(ii)=sum(isq(:));
end

