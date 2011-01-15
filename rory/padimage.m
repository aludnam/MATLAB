function imout=padimage(im, sizeimout)

sizeim = size(im);

sizedifhalf = floor((sizeimout-sizeim)/2);

imout = padarray(im, sizedifhalf);

sizedif = sizeimout-size(imout);

if sizedif(1)>0
    imout = [imout; zeros(1,size(imout,2))];
end
if sizedif(2)>0
    imout = [imout, zeros(size(imout,1),1)];
end