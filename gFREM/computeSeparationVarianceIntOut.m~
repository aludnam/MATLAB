
function [vardintout, Iintout]=computeSeparationVarianceIntOut(x,l1,l2,sig,int_vec, pixelizeversion, offset)
% [vardintout, Iintout]=computeSeparationVarianceIntOut(x,l1,l2,sig,int_vec, pixelizeversion)
% Computes variance (as a inverse of the Fisher Information) of two points
% separated by a distance d=c1-c2. Integrating out the intensity. 

if ~exist('offset','var')
    offset=0; % background
end

f1=makeGauss(x,l1,sig(1));                  % creates PSF (gauss approx)
for ind_dist=1:length(l2)                   % distance    
    f2=makeGauss(x,l2(ind_dist),sig(2));    % creates PSF (gauss approx) shifted to l2       
    nphot = int_vec(1);    
    Iintout(:,:,ind_dist)=numericalMeanEstimation(x,f1,f2, offset, nphot);
    vardintout(ind_dist)=[1,-1]/Iintout(:,:,ind_dist)*[1,-1]';
end
q=0;
