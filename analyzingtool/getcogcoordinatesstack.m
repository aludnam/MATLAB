function [coord, inbin]=getcogcoordinatesstack(in,thresholdmax)
% [coord, inbin]=getcoordinates(in,thresholdmax)
% gets coordinates of center of gravity of hte clusters from the (CSSTORM) stack of images images. 
% in - input image (grey value)
% thresholdmax - threshold wrt maximum in each frame to be considered as a
% valid local maximum (default: thresholdmax = .1)
% coord(:,[1 2]) -  coordinates of hte cog
% coord(:,3) - slice of hte corresponding cluster
% coord(:,4) - value of the cluster maximum
% inbin - binary threholded input

if ~exist('thresholdmax','var')
    thresholdmax=0.1; 
end

s=size(in); 
invec=reshape(in,s(1)*s(2),s(3)); 
invecn=normcMax(invec); 
inbinvec=invecn>thresholdmax;
inbin=reshape(inbinvec,s(1),s(2),s(3)); 

coord=[];
for ii=1:s(3); % unforunately i have to do this rahter then calling the gravity as i can't get the right connectivity for the stack of images...
    data= measure(squeeze(inbin(:,:,ii)),squeeze(in(:,:,ii)),{'gravity','maxval'},[],2);
    d=data.Gravity';
    d(:,3)=ii;
    dm=[d,data.maxVal'];
    coord=[coord;dm];    
end