function simNoise = generateNoise(nPoints, sizeVec)
% simNoise = generateNoise(nPoints, sizeVec)
% 
% generates Poisson noise with given number of points (nPoints) and given
% image size (sizeVec = [xlim1, xlim2, ylim1, ylim2]).

simNoise(:,1) = sizeVec(1) + rand(nPoints,1)*(sizeVec(2) - sizeVec(1));
simNoise(:,2) = sizeVec(3) + rand(nPoints,1)*(sizeVec(4) - sizeVec(3));