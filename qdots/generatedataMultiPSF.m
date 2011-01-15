function [dpixc, dpixc_ind, blinkmat] = generatedataMultiPSF(size_vec, x_vec,y_vec,intensity_vec,psf_tmp,offset,Nt, probtrans)
% generate points

if size(intensity_vec,2)>size(intensity_vec,1); intensity_vec = intensity_vec'; end
if size(x_vec,2)>size(x_vec,1); x_vec = x_vec'; end
if size(y_vec,2)>size(y_vec,1); y_vec = y_vec'; end

nx = size_vec(2)-size_vec(1);
ny = size_vec(4)-size_vec(3);
N = length(x_vec); %number of points
center = [nx ny]/2;
center_round = round(center);
 
% 
% intensity_mat = repmat(intensity_vec, 1,Nt);
% changemat = rand(N,Nt)<probtrans;
% statemat = mod(cumsum(changemat,2),2);
% initvec = rand(N,1)>0.5; %initial state of the blinkmat
% ivt = ~(initvec == statemat(:,1));
% statematinit = mod(statemat+repmat(ivt,1,Nt),2);
% 
% % blinkmat_equal = rand(N, Nt);
% % blinkmat = blinkmat_equal .* intensity_mat; %different intensities...
% blinkmat = statematinit .* intensity_mat; %different intensities...
blinkmat = blinkmat_markov(N,Nt, intensity_vec, probtrans);
center_im = pixelize(center_round, 1, size_vec, nx, ny, [],0);
dpixc_ind = newimar(N);

shift_vec = [x_vec-center(1), y_vec-center(2)];
for ii=1:N
    sp = size(psf_tmp{ii});
    psf_tmp2 = double(psf_tmp{ii}/sum(psf_tmp{ii}(:))); %normalize to 1
    psf_tmp3 = padimage(psf_tmp2,fliplr([nx, ny])+2*sp);
    psf = psf_tmp3;
    psf_im = clip(dip_image(conv2(center_im,psf,'same')), 0, Inf);
    dpixc_ind{ii} = clip(shift(psf_im, shift_vec(ii,:)), 0, Inf);
end

dpixc_nonoise = array2im(dpixc_ind'*blinkmat + offset);
dpixc_dip = noise(dpixc_nonoise,'poisson');
dpixc = double(dpixc_dip);

