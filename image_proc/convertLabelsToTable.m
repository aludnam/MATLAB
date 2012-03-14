function t=convertLabelsToTable(L,n)
%convertLabelsToTable Converts labeled stack of images to the table of cluster-sizes.
%   t=convertLabelsToTable(L,n)
%   L is a labelled image, n number of clustrers - output of the function
%   [L,n] = bwlabelStack(bw,ncon);
%   t is the table (#frames X #clusters) containing sizes of the individual
%   clusters in each frame.
sizeZ=size(L,3); 
maxClust = max(n); 
t=zeros(sizeZ,maxClust);
for frame=1:sizeZ
    for cluster = 1:maxClust
        clustTmp=L(:,:,frame)==cluster;
        t(frame,cluster)=sum(clustTmp(:));
    end
end
