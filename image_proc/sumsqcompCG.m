function sumsq = sumsqcompCG(dpixc, res, htrace1, htrace2, peval)
ROIx=peval.ROIx(1):peval.ROIx(2);
ROIy=peval.ROIy(1):peval.ROIy(2);
ROIz=peval.ROIz(1):peval.ROIz(2);
dpixcs = dpixc(ROIx,ROIy,ROIz);
dveccs = reshape(dpixcs, peval.numpix, peval.nt);
htr =exp(res.htrace);
for jj=1:size(htrace1,1) %avg
    for ii=1:size(htrace1,3) %iter
        
        htmp=[squeeze(exp(htrace1(jj,:,ii))); squeeze(exp(htrace2(jj,:,ii))); res.h(3,:)];
        itmp=res.w*htmp;
        
        isq=(dveccs-itmp).^2;
        sumsq(ii,jj)=sum(isq(:));
    end
end
