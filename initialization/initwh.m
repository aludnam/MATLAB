function [winit, hinit] = initwh(winit_pix, hinit, peval, verbose)
% [winit, hinit] = initwh(winit_pix, hinit, peval, verbose)
if ~exist('verbose','var')
    verbose=1;
end

if isempty(winit_pix) %random initialization
    if verbose
        fprintf('Uniform random initialization of W.\n')
    end
    winit = rand(peval.nx*peval.ny, peval.ncomp);
else
    switch class (winit_pix)
        case 'dip_image_array'
            winit_pix = double(array2im(winit_pix));
        case 'dip_image'
            winit_pix = double(winit_pix);
    end
    winit = reshape(winit_pix, peval.nx*peval.ny, peval.ncomp);
end

if isempty(hinit) %random initialization
    if verbose
        fprintf('Uniform random initialization of H.\n')
    end
    hinit = rand(peval.ncomp, peval.nt);
end
    
w_bg = ones(peval.nx*peval.ny, 1); %background component
winit_bg = [winit, w_bg];
sumw = sum(winit_bg,1);
winit = winit_bg./repmat(sumw, peval.nx*peval.ny, 1); %normalized to 1

h_bg = peval.bg*peval.nx*peval.ny*ones(1, peval.nt); %background component

fm = winit(:, 1:peval.ncomp)*hinit;
f=mean(fm(:));

hinit = [(peval.meandata-peval.bg)/f * hinit; h_bg]; %to get same mean as data

if verbose
    fprintf('One component [%g] added as a background...\n', peval.ncomp+1)
end