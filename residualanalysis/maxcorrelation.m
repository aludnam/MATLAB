function maxcorr=maxcorrelation(res, dpixc, peval, linkagemethod, correlationmethod)
% maxcorr=maxcorrelation(res, dpixc, peval, linkagemethod, correlationmethod)
% Computes maximum correlation of the normalized residuals
% linkagemethod: method for liknage function [>help linkage], default =
% 'average'
% correlationmethod: method for computing distance from correlation matrix
% [>help correlation2distance], default = 'pos' 

% Default values:
default_linkagemethod = 'average';
default_correlationmethod = 'pos';

if ~exist('linkagemethod','var')
    linkagemethod = default_linkagemethod;
    fprintf('Default linkage method: %s\n', linkagemethod)
end
if ~exist('correlationmethod','var')
    correlationmethod = default_correlationmethod;
    fprintf('Default correlation method: %s\n', correlationmethod)
end

% Computes residuals:
[resid, resid_norm, Wxk,Hkt,centers,Vxtpix, Vxtpixbg] = computeresid(res, dpixc, peval);
data=resid_norm;
sized = size(data);
dveccr= reshape(data,sized(1)*sized(2), sized(3));

% Computes correlation coeeficients:
ccd = (corrcoef(dveccr'));

% Transforms coorelations to 'distance matrix' 
ccds = correlation2distance(ccd, correlationmethod);

% Hierarcical cluster tree. (Dendrogram by dendrogram(Z))   
Z = linkage(ccds,linkagemethod);                
maxcorr=1-min(Z(:,3));