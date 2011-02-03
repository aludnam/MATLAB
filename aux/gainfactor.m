function c = gainfactor(EMGain)
% c = gainfactor(EMGain)
% computes theoretical conversion rate for photon numbers
% number of photons = c * intensity 
% EMGain =130;
b = 4*exp(0.02223*EMGain);
c = 5.8/(0.92*b);