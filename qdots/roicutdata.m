dpixc = dpixc(ROIx,ROIy,ROIz);
if exist('dpixc_ind','var');
ditmp = double(array2im(dpixc_ind));
ditmp = ditmp(ROIx,ROIy,:);
dpixc_ind = im2array(dip_image(ditmp));
end