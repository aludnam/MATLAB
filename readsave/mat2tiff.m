function mat2tiff(image,namedir, namefilebase)
% mat2tiff(image,namedir, namefilebase)
% Saves image (2D image or 3D stack) as a series of 16-bit tiff images into direcotry namedir.
% 
% Example:  mat2tiff(dpixc,'images','dpixc')

s=size(image);
mkdir(namedir)
if ~strcmpi(class(image),'uint16')
    image=uint16(image);
    warning('Converting to uint16 type!')
end

if numel(s)==3
    for ii=1:s(3)
        namefile = [namefilebase '_' sprintf('%0.4d',ii)];
        imwrite(image(:,:,ii),[namedir '/' namefile '.tiff'] ,'TIFF')
    end
elseif numel(s)==2
    imwrite(image,[namedir '/' namefilebase '.tiff'] ,'TIFF')
else
    errror('Wrong dimensionality of the input figure!')
end
