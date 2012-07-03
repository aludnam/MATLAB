function out=fliplrstack(in)
% out=fliplrstack(in)
% Uses fliplr on each frame of the stack in. 

isdip = 0; 
if strcmpi(class(in),'dip_image')
    isdip=1; 
    fprintf('Converting the input image temporarly to double\n')
    in = double(in); 
    
end

sz=size(in,3);
out = zeros(size(in));
for ii=1:sz
    out(:,:,ii)=fliplr(in(:,:,ii)); 
end

if isdip 
    out=dip_image(out);
    fprintf('Converting back to dip_image.\n')
end

