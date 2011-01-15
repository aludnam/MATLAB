function [xK, K, KNMax, KNMin] = ...
    kfunction_main (dataXY_all, xlim1, xlim2, ylim1, ylim2, Klim, nSteps, envelopes, Nsimul, filename)

% KFUNCTION_MAIN computes Ripley's K function for soelected ROI in the data
% [xK, K, KNMax, KNMin] = ...
%     kfunction_main (dataXY_all, xlim1, xlim2, ylim1, ylim2, Klim, nSteps,
%     envelopes, Nsimul, filename)
% dataXY_all - input data, #points-by-2 matrix, each row correspond to xy
% coordinate of the datapoint
% xlim1, xlim2, ylim1, ylim2 - selected ROI
% Klim - limits of the K function ([Klim_low, Klim_high]) should be less
% then half of the smaler size of the ROI
% nSteps - number of steps for computation of K function
% envelopes - if set to 1 computes simulated envelopes for the same number
% of datapoints
% Nsimul -  number of simulations to compute envelops
% filename - name of the file where data are written (in current director)
% If ommited no file is created.
p.xlim1=xlim1; p.xlim2=xlim2;  p.ylim1=ylim1; p.ylim2=ylim2; 
p.Klim=Klim; p.nSteps=nSteps; p.envelopes=envelopes;

dataXY = ROIdata(dataXY_all, p.xlim1, p.xlim2, p.ylim1, p.ylim2, 0);
p.nPoints = length(dataXY);

p.xStep  = (p.Klim(2)-p.Klim(1))/p.nSteps;
xK = (0:p.xStep: round(p.xStep*p.nSteps))';

box = [p.xlim1, p.xlim2, p.ylim1, p.ylim2];
fprintf('Computing Ripley''s K-function for selected ROI... \n');
K = kfunction(dataXY, xK, box, 1);


if p.envelopes
    p.Nsimul=Nsimul;
    fprintf('Computing simulation envelopes... \n');
    for ii=1:p.Nsimul
        simNoise = generateNoise(p.nPoints,box);
        KN = kfunction(simNoise,xK, box, 1);
        if ii==1 %first round
            KNMax = KN;
            KNMin = KN;
        else
            KNMax = max(KNMax, KN);
            KNMin = min(KNMin, KN);
        end
    end
    if nargin>9
        writedata(xK, computeL(xK, KNMax), p,  [filename '_envelope_Max'], ['L-function - Maximum of ' num2str(p.Nsimul) ' simulations of Poisson process']);
        writedata(xK, computeL(xK, KNMin), p,  [filename '_envelope_Min'], ['L-function - Minimum of ' num2str(p.Nsimul) ' simulations of Poisson process']);
    end
else
    KNMax = [];
    KNMin = [];
end

if nargin>9
    writedata(xK, computeL(xK, K), p,  filename, 'L-function')
    save (filename)

%     writedata(xK, K, KNMax, KNMin, p,  filename, 'Ripley''s K-function')
end
