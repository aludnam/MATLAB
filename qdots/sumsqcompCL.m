function sumsq = sumsqcompCL(dpixc, res, htrace1, htrace2, peval)
ROIx=peval.ROIx(1):peval.ROIx(2);
ROIy=peval.ROIy(1):peval.ROIy(2);
ROIz=peval.ROIz(1):peval.ROIz(2);
dpixcs = dpixc(ROIx,ROIy,ROIz); 
dveccs = reshape(dpixcs, peval.numpix, peval.nt);
htr =res.htrace;
for jj=1:size(htrace1,1)    %avgs
    for ii=1:size(htrace1,3) %iter
%         itmp=res.w*squeeze(res.htrace(ii,:,:));
        itmp=res.w*[squeeze(htrace1(jj,:,ii)); squeeze(htrace2(jj,:,ii)); squeeze(res.h(3,:))];
        isq=(dveccs-itmp).^2;
        sumsq(ii,jj)=sum(isq(:));
    end
end
