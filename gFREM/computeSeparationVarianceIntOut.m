
function [vardintout, Iintout]=computeSeparationVarianceIntOut(x,xhires,l1,l2,sig,int_vec, pixelizeversion, offset)
% [vardintout, Iintout]=computeSeparationVarianceIntOut(x,l1,l2,sig,int_vec, pixelizeversion)
% Computes variance (as a inverse of the Fisher Information) of two points
% separated by a distance d=c1-c2. Integrating out the intensity. 

if ~exist('offset','var')
    offset=0; % background
end
oversampleFactor=round((x(1,2,1)-x(1,1,1))/(xhires(1,2,1)-xhires(1,1,1)));
if ndims(x) == 2 %1D vector
    f1=normcSum(makeGauss(x,l1,sig(1)));                  % creates PSF (gauss approx -> !!! Different to simulationtools/makegauss.m !!!)
else 
    f1_2D=makeGauss2D(xhires,l1,sig(1));
    f1 = normcSum(reshape(f1_2D,1,numel(f1_2D)));          % making 1D vector by concatenating 2D array
end
    
for ind_dist=1:length(l2)                       % distance    
    if ndims(x)==2 %1D vector
        f2=nomrcSum(makeGauss(x,l2(ind_dist),sig(2)));    % creates PSF (gauss approx) shifted to l2
    else
        f2_2D=makeGauss2D(xhires,l2(ind_dist),sig(2));
        f2 = normcSum(reshape(f2_2D,1,numel(f2_2D)));     % making 1D vector by concatenating 2D array
    end
    nphot = int_vec(1);    
    Iintout(:,:,ind_dist)=numericalMeanEstimation(xhires,f1,f2, offset, nphot,oversampleFactor);
    vardintout(ind_dist)=[1,-1]/Iintout(:,:,ind_dist)*[1,-1]';
end
q=0;
