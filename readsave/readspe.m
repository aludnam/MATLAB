function imout=readspe(filename)
% imout=readspe(filename)
% Reads .spe file from the CSSTORM evaluation
fid = fopen(filename,'r','l');
header = fread(fid,2050,'*uint16');
imgs = fread(fid, inf, 'float32');
fclose(fid);
z_dim = double(header(724));
x_dim = double(header(22));
im=reshape(imgs,x_dim,x_dim,z_dim);
imnan=isnan(im);
imn=zeros(size(im));
imn(~imnan)=im(~imnan);
imout=rot90stack(fliplrstack(imn));
