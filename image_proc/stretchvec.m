function out=stretchvec(in,minimum,maximum)
% out=stretchvec(in,minimum,maximum)
% Stretches vector in between minimum and maximum value. 
% Defaults      minimum = 0; 
%               maximum = 1;   
if ~exist('minimum','var'); minimum = 0; end
if ~exist('maximum','var'); maximum = 1; end
minin=min(in); 
maxin=max(in); 

normin=(in-minin)/(maxin-minin); % between 0 and 1
out = normin*(maximum-minimum)+minimum; 