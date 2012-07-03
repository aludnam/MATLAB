function out=maximastack(A,connectivity)
% out=maximastack(A,connectivity)
% Gives  binary image of local maxima of each frame in the stack A. Type
% dip_maxima for more information about connectivity (either 1 or 2). 


if ~exist('connectivity','var')
    connectivity=2; % 8 connected neigbours in 2d... Set to 1 to get 4 connected neigbours only. 
end
out=zeros(size(A)); 
for ii=1:size(A,3)    
    out(:,:,ii) = dip_maxima(A(:,:,ii),[],connectivity,1);
%     m{ii}=sort(A(maxpix),'descend'); 
end