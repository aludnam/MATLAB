function X = sampleimages( winsize, samples )
% sampleimages - gathers image patches from Olshausens images.
%

% Start by loading the images
fprintf('Loading images...\n');
load ../data/IMAGES;
num_images=size(IMAGES,2);
image_size=sqrt(size(IMAGES,1));

% This will hold the patches
X=zeros(winsize^2,samples);
totalsamples = 0;

for i=1:num_images
  
  % Choose an image for this batch
  this_image=reshape(IMAGES(:,i),image_size,image_size)';
  BUFF=4;

  % Determine how many patches to take
  getsample = floor(samples/num_images);
  if i==num_images, getsample = samples-totalsamples; end
  
  % Extract patches at random from this image to make data vector X
  for j=1:getsample
    r=BUFF+ceil((image_size-winsize-2*BUFF)*rand);
    c=BUFF+ceil((image_size-winsize-2*BUFF)*rand);
    totalsamples = totalsamples + 1;
    X(:,totalsamples) = ...
	reshape( this_image(r:r+winsize-1,c:c+winsize-1),winsize^2,1);
  end
  
end  
