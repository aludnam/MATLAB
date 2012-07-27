function edgedist=maxtoedgedist(im)
% edgedist=maxtoedgedist(im)
% Computes the distance of hte global maximum from teh edge. 
sx=size(im,1); 
sy=size(im,2); 
sz=size(im,3);

edgedist=inf(sz,1); 
for ii=1:sz
    imframe=im(:,:,ii); 
    mim=max(max(imframe)); 
    [maxx,maxy]=find(imframe==mim,1);
    edgedist(ii)=min([maxx-1,sx-maxx, maxy-1, sy-maxy]); 
end