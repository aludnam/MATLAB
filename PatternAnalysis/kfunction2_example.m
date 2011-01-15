function kfunction2_example
% example of K-function2 analysis of the data


%load synthetic clustered data(simData) (molecules)
load simData.mat
% or random position of simData (uncomment following line...)
% simData = generateNoise(1000, [0, 60, 0, 60]);

%load centers of clusters (vesicles)
load meanClust.mat
% or random position of meanClust (uncomment following line...)
% meanClust = generateNoise(20, [0 60 0 60]);

% plots simData
plotData(simData, [0 60 0 60])

% plots meanClust on top of simData in red color in size 20
colorData = 'r'; sizeData = 20;
hold on
plotData(meanClust, [0 60 0 60] ,colorData, sizeData)
hold off
xlim1 = 0; xlim2 = 60; ylim1 = 0; ylim2 = 60;
% computes K function in for values between Klim with nSteps steps for ROI
% and simulation envelopes for
Klim = [0 10]; % range of the K function
nSteps = 100; % step of the position of the K function esstimation
envelopes_method = 1;  %to compute simulatyion envelopes from random data1 and data2
Nsimul = 3; % number of simulations for envelope coputation
filename = 'Kfunction2_example'; % name of the file where K function values are stored
 
[xK, K, KNMax, KNMin] = ...
    kfunction2_main(meanClust, simData,xlim1, xlim2, ylim1, ylim2, Klim, nSteps,...
    envelopes_method, Nsimul, filename); 

% plot K function with theoretical line for Poisson process K_p = pi*xK^2:
plot (xK, K, xK, pi*xK.^2, '--k')

% plot simulation envelopes:
hold on 
plot (xK, KNMax, '-.r', xK, KNMin, '-.r')

% compute L function 
L = computeL(xK, K);

% plot L function with computed envelopes
figure
plot(xK, L, xK, computeL(xK, KNMax), '-.r', xK, computeL(xK, KNMin), '-.r')
hold on
plot([xK(1), xK(end)], [0 0], '--k')