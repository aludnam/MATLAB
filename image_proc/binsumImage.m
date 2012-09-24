function out = binsumImage(in, bin)
% Rebins image by summing bin(1)xbin(2) pixels together along x and y
% dimension, respectively. Starts at upper left corner. Rebinned image will
% have floor(size(in)./bin) pixels.
%
% out = binsum(in, bin)
di=0;
if strcmp(class(in),'dip_image')
    di=1;
    in=double(in);
end

npix=ceil(size(in)./bin);
out=zeros(npix); 
for xi=1:npix(1)
    startx=(xi-1)*bin(1)+1;    
    endx=min(startx+bin(1)-1,size(in,1)); 
    for yi=1:npix(2)
        starty=(yi-1)*bin(2)+1;                
        endy=min(starty+bin(2)-1,size(in,2));
        binsum=sum(sum(in(startx:endx,starty:endy)));                 
        out(xi,yi)=binsum;
    end    
end

if di
    out=dip_image(out); 
end
