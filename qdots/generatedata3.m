function [dpixc, dveccr, dpixc_ind, blinkmat, N] = generatedata3(sizevec, psf, separ, maxphotvec, offset, Nt)
% generate 2 points with different shift, diffrent PSFs

nx = sizevec(2)-sizevec(1);
ny = sizevec(4)-sizevec(3);
N = 2;

if size(maxphotvec,2)>size(maxphotvec,1); maxphotvec = maxphotvec'; end
maxphotmat = repmat(maxphotvec, 1,Nt);

blinkmat_equal = rand(N, Nt);
blinkmat = blinkmat_equal .* maxphotmat; %different intensities...
center = round([nx ny]/2);
centerim = pixelize(center, 1, sizevec, nx, ny, [],0);
dpixc_ind = newimar(2);
dpixc_ind{1} = clip(dip_image(conv2(centerim,psf{1},'same')), 0, Inf);
dpixc_ind_tmp = clip(dip_image(conv2(centerim,psf{2},'same')), 0, Inf);

dpixc_ind{2} = clip(shift(dpixc_ind_tmp, separ), 0, Inf);

dpixc_nonoise = array2im(dpixc_ind'*blinkmat + offset);

dpixc_dip = noise(dpixc_nonoise,'poisson');
dpixc = double(dpixc_dip);
dveccr = double(squeeze(reshape(dpixc, nx*ny, 1, Nt))); % vectors of resized images

