function pow=makePowerFromLabel(labelIn, defPow)
% pow=makePowerFromLabel(labelIn, defPow)

pow = ones(size(labelIn)); 
pow(labelIn==1)=defPow;
pow(labelIn==2)=defPow;
