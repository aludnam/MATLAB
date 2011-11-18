function [dpixc, dpixc_ind, blinkmat, psf_z_index] = generatedataMultiBlinkmat(size_vec, x_vec,y_vec,psf_tmp,offset,Nt, blinkmat,noisetype)
% [dpixc, dpixc_ind, blinkmat, psf_z_index] = generatedataMultiBlinkmat(size_vec,
% x_vec,y_vec,psf_tmp,offset,Nt, blinkmat,noisetype)
% generate points
if ~exist('noisetype','var')
    noisetype='Poisson';
end

if size(x_vec,2)>size(x_vec,1); x_vec = x_vec'; end
if size(y_vec,2)>size(y_vec,1); y_vec = y_vec'; end

nx = size_vec(2)-size_vec(1);
ny = size_vec(4)-size_vec(3);
N = length(x_vec); %number of sources
center = [nx ny]/2;
center_round = round(center);
psf_z_size = size(psf_tmp,3); 
center_im = zeros(nx, ny, psf_z_size);
center_im(center_round(1), center_round(2), :)=1; 
dpixc_ind = newimar(N);

psf_double = double(psf_tmp); % to unify indexing (starting from 1)
psf_middle = psf_double(:,:,ceil(psf_z_size/2));
psf = (psf_tmp/sum(psf_tmp(:))); %normalize so that in-focus (middle) plane sums to 1

psf_im = newim(nx, ny, psf_z_size); % this will be PSF centerd and at the required size and in dip_image format
for ii=1:psf_z_size
    p = clip(dip_image(conv2(center_im(:,:,ii),psf(:,:,ii),'same')), 0, Inf); %makes the PSF the same size and places in the round(center)
    psf_im(:,:,ii-1) = clip(shift(p, center - center_round+0.5)); %shift to the proper center
end
shift_vec = [x_vec-center(1), y_vec-center(2)]+.5;

psf_z_index=randi(size(psf_im,3),1,N)-1; % indices of the z-plane for individual PSF used for simulation

for ii=1:N
    if ndims(psf_im)==2
        psf_use = squeeze(psf_im(:,:,psf_z_index));
    elseif ndims(psf_im)==3
        psf_use = squeeze(psf_im(:,:,psf_z_index(ii)));
    end
    dpixc_ind{ii} = clip(shift(psf_use, shift_vec(ii,:)), 0, Inf);
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


