function winit_pix = resetwinit(winit_pix, reset_vec, peval);
for ii=1:length(reset_vec)
    winit_pix(:,:,reset_vec(ii))=rand(peval.nx,peval.ny);
end