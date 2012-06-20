function C=convstack(A,B,shape)
% C=convstack(A,B,shape) makes convolution of the each slice of the stack of images A with B.
% It is calling conv2 funciton (help conv2)
%   C = CONV2(..., SHAPE) returns a subsection of the 2-D
%     convolution with size specified by SHAPE:
%       'full'  - (default) returns the full 2-D convolution,
%       'same'  - returns the central part of the convolution
%                 that is the same size as A.
%       'valid' - returns only those parts of the convolution
%                 that are computed without the zero-padded
%                 edges. size(C) = [ma-mb+1,na-nb+1] when
%                 all(size(A) >= size(B)), otherwise C is 
%                 an empty matrix [].

if ~exist('shape','var')
    shape = 'same';
end

C=zeros(size(A));
for ii=1:size(A,3)
   C(:,:,ii)=conv2(A(:,:,ii),B,shape);
end
    
    

