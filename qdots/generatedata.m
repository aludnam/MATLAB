function [dpixc, dveccr, N] = generatedata(d, sizevec, psf, maxphot, offset, Nt, rs)

nx = sizevec(2)-sizevec(1);
ny = sizevec(4)-sizevec(3);
N = size(d,1);

blinkmat = rand(N, Nt);

for ii=1 : N
    dpix_ind(:,:,ii) = pixelize(d(ii,:), 1, sizevec, nx, ny, [],0);
    dpixc_ind(:,:,ii) = conv2(dpix_ind(:,:,ii),psf,'same'); %convolution of individual points    
end

dvecc_ind = squeeze(reshape(dpixc_ind, nx*ny, 1, N)); % each indiv image to vector

dvecc_nonoise = (dvecc_ind*blinkmat); % set of Nt vectors
dpixc_nonoise = imresize(reshape(dvecc_nonoise, nx, ny, Nt), rs);
dpixc_nonoise_dip = dip_image(dpixc_nonoise);

maxdvec = max(dpixc_nonoise(:));
offset_abs = offset/(1-offset)*maxdvec;
dpixc_nonoise_dip = (dpixc_nonoise_dip + offset_abs)/(maxdvec+offset_abs)*maxphot;
dpixc_dip = noise(dpixc_nonoise_dip,'poisson');
dpixc = double(dpixc_dip);
dveccr = double(squeeze(reshape(dpixc, nx*ny*rs^2, 1, Nt))); % vectors of resized images

