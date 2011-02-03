function [dpixc, dpixc_ind, blinkmat] = generatedataMultiBlinkmat(size_vec, x_vec,y_vec,intensity_vec,psf_tmp,offset,Nt, blinkmat,noisetype)
% generate points
if ~exist('noisetype','var')
    noisetype='Poisson';
end
if size(intensity_vec,2)>size(intensity_vec,1); intensity_vec = intensity_vec'; end
if size(x_vec,2)>size(x_vec,1); x_vec = x_vec'; end
if size(y_vec,2)>size(y_vec,1); y_vec = y_vec'; end

nx = size_vec(2)-size_vec(1);
ny = size_vec(4)-size_vec(3);
N = length(x_vec); %number of points
center = [nx ny]/2;
center_round = round(center);
 
sp = size(psf_tmp);
psf_tmp2 = double(psf_tmp/sum(psf_tmp(:))); %normalize to 1
% psf_tmp3 = padimage(psf_tmp2,fliplr([nx, ny])+2*sp);
psf = psf_tmp2;

%intensity_mat = repmat(intensity_vec, 1,Nt);
%blinkmat_equal = rand(N, Nt);
%blinkmat = blinkmat_equal .* intensity_mat; %different intensities...
center_im = pixelize(center_round, 1, size_vec, nx, ny, [],0);
dpixc_ind = newimar(N);
psf_im = clip(dip_image(conv2(center_im,psf,'same')), 0, Inf); %makes the PSF the same size and places in the round(center)

psf_im = clip(shift(psf_im, center - center_round+0.5));

shift_vec = [x_vec-center(1), y_vec-center(2)]+.5;

for ii=1:N
    dpixc_ind{ii} = clip(shift(psf_im, shift_vec(ii,:)), 0, Inf);
end

dpixc_nonoise = array2im(dpixc_ind'*blinkmat + offset);
if ~isempty(noisetype)
    fprintf('Adding %s noise.\n', noisetype)
    dpixc_dip = noise(dpixc_nonoise,noisetype);
else
    fprintf('No noise added!\n');
    dpixc_dip = double(dpixc_nonoise);
end
dpixc = double(dpixc_dip);
% fprintf('No noise!');
% dpixc = double(dpixc_nonoise);


