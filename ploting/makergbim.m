function h=makergbim(im1, im2, channel)
% h=makergbim(im1, im2, channel)
% Creates RGB image where im1 is as a base image (white) and im2 is
% overlayed over channel ([1 0 0] for red, [0 1 0] for green....).
% im1 and im2 must be the same size.
% h is a handle to the image


if ~exist('channel', 'var')
    channel = [1 0 0];
end

m1=max(im1(:));
m2=max(im2(:));
rm=m1/m2;

jc=joinchannels('rgb', im1+channel(1)*rm*im2, im1+channel(2)*rm*im2, im1+channel(3)*rm*im2);
h=dipshow(jc);