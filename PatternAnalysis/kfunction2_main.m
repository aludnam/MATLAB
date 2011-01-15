function [xK, K, KNMax, KNMin] = ...
    kfunction2_main (data1_all, data2_all, xlim1, xlim2, ylim1, ylim2, Klim, nSteps, envelopes_method, Nsimul, filename)

% KFUNCTION2_MAIN computes Ripley's K function for soelected ROI in the data
% [xK, K, KNMax, KNMin] = ...
%     kfunction2_main (data1_all, data2_all, xlim1, xlim2, ylim1, ylim2, Klim, nSteps,
%     envelopes, Nsimul, filename)
% data1, data2 - input data, #points-by-2 matrix, each row correspond to xy
% coordinate of the datapoint
% xlim1, xlim2, ylim1, ylim2 - selected ROI
% Klim - limits of the K function ([Klim_low, Klim_high]) should be less
% then half of the smaler size of the ROI
% nSteps - number of steps for computation of K function
% envelopes_method : computes simulation envelopes
%     envelopes_method = 0 -> no envelopes computed
%     envelopes_method = 1 -> simulates random data1 but data2 leaves as original
%     envelopes_method = 2 -> simulates random data2 but data1 leaves as original
%     envelopes_method = 3 -> simulates random data1 and data2 (same number of points)
% Nsimul -  number of simulations to compute envelops
% filename - name of the file where data are written (in current director)
% If ommited no file is created.

p.xlim1=xlim1; p.xlim2=xlim2;  p.ylim1=ylim1; p.ylim2=ylim2;
p.Klim=Klim; p.nSteps=nSteps; p.envelopes_method=envelopes_method; p.Nsimul=Nsimul;

data1 = ROIdata(data1_all, p.xlim1, p.xlim2, p.ylim1, p.ylim2, 0);
data2 = ROIdata(data2_all, p.xlim1, p.xlim2, p.ylim1, p.ylim2, 0);

p.nPoints1 = length(data1);
p.nPoints2 = length(data2);

p.xStep  = (p.Klim(2)-p.Klim(1))/p.nSteps;
xK = (0:p.xStep: round(p.xStep*p.nSteps))';

box = [p.xlim1, p.xlim2, p.ylim1, p.ylim2];
fprintf('Computing Ripley''s K-function for selected ROI... \n');
K = kfunction2(data1, data2, xK, box, 1);

if p.envelopes_method > 0
    fprintf('Computing simulation envelopes... \n')
    switch envelopes_method
        case 1
            fprintf('Simulating random data1.\n');
        case 2
            fprintf('Simulating random data2.\n');
        case 3
            fprintf('Simulating random data1 and data2.\n');
    end
    for ii=1:p.Nsimul
        switch envelopes_method
            case 1
                data1 = generateNoise(p.nPoints1,box);
            case 2
                data2 = generateNoise(p.nPoints2,box);
            case 3
                data1 = generateNoise(p.nPoints1,box);
                data2 = generateNoise(p.nPoints2,box);
        end
        KN = kfunction2(data1, data2,xK, box, 1);
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
