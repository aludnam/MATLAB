function V = natdata
% natdata - read natural image data preprocessed into ON/OFF channels
%

% Data parameters
wsize = 12;
samples = 10000;
    
% Sample image patches from pre-whitened Olshausen's images
Z = sampleimages( wsize, samples ); 
   
% Separate into ON and OFF channels.
Y = Z; Y(Y<0)=0; Z(Z>0)=0; Z=-Z;
V = [Y; Z];

% Additionally, this is required to avoid having any exact zeros:
% (divergence objective cannot handle them!)
V = max(V,1e-4);

% Normalize each channel to unit mean squared activation
V = V./(sqrt(mean(V.^2,2))*ones(1,size(V,2)));
    

