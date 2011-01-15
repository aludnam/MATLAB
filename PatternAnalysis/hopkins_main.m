function [xH, H, HNMax, HNMin] = hopkins_main (dataXY_all, xlim1, xlim2, ylim1, ylim2, m, N, Nbins, envelopes, Nsimul, filename)
% HOPKINS_MAIN computes hopkins statistics
% [xH, H, HNMax, HNMin] = hopkins_main (dataXY_all, xlim1, xlim2, ylim1, ylim2, m, N, Nbins, envelopes, Nsimul, filename)
% dataXY_all - whole data
% [xlim1, xlim2, ylim1, ylim2] - selected ROI
% m - number of random events and points for Hopkins computation 
% N - number of iterations for histogram computation 
% Nbins - numbero of bins in histogram 
% envelopes - if set to 1 maximum and minimum envelopes for Poisson random 
% process are computed
% Nsimul - number of simulations for envelopes computation 
% filename - name of the txt file where results are written 

Hopp.xlim1=xlim1; Hopp.xlim2=xlim2; Hopp.ylim1=ylim1; Hopp.ylim2=ylim2;
Hopp.Nbins = Nbins; Hopp.m = m; Hopp.N = N; 
if envelopes, Hopp.Nsimul = Nsimul; end

figure
dataXY = ROIdata(dataXY_all, Hopp.xlim1, Hopp.xlim2, Hopp.ylim1, Hopp.ylim2, 1); 
Hopp.Npoints = length(dataXY);

fprintf('Computing Hopkins statistics for selected ROI... \n');
h = waitbar(0,'Computing Hopkins...');
for ii=1:Hopp.N
    HopD(ii) = hopkins(dataXY, Hopp.xlim1, Hopp.xlim2, Hopp.ylim1, Hopp.ylim2, Hopp.m);
    waitbar(ii/Hopp.N,h)
end
close (h)

[H, xH] = hist(HopD, Nbins);
H = H/trapz(xH, H); % normalization

if envelopes
    [HNMin, HNMax, xH] = ...
        hopkinsEnv(Hopp.Npoints, Hopp.xlim1, Hopp.xlim2, Hopp.ylim1, Hopp.ylim2, Hopp.m,Hopp.N, xH, Hopp.Nsimul);
    if nargin>10
        writedata(xH, HNMax, Hopp,  [filename '_envelope_Max'], ['Hopkins statistics - Maximum of ' num2str(Hopp.Nsimul) ' simulations of Poisson process']);
        writedata(xH, HNMin, Hopp,  [filename '_envelope_Min'], ['Hopkins statistics - Mainimum of ' num2str(Hopp.Nsimul) ' simulations of Poisson process'])
    end
else
    HNMin = [];
    HNMax = [];
end

if nargin>10
    writedata(xH, H, Hopp,  filename, 'Hopkins statistics')
    save (filename)
end 


function [HhistMin, HhistMax, xHhist] = hopkinsEnv(Npoints, xlim1, xlim2, ylim1, ylim2, m,N, bins, Nsimul)
% function [HhistMin, HhistMax, xHhist] = hopkinsEnv(Npoints, xlim1, xlim2,
% ylim1, ylim2, m,N, bins, Nsimul)

h = waitbar(0,'Computing Hopkins Envelope...');
for ii=1 : Nsimul
    waitbar(ii/Nsimul,h)
    simNoise = generateNoise(Npoints, [xlim1,xlim2,ylim1,ylim2]);
    for jj=1 : N
        H(jj) = hopkins (simNoise,xlim1, xlim2, ylim1, ylim2, m);
    end
    [Hhist, xHhist] = hist(H, bins);
    Hhist = Hhist/trapz(xHhist, Hhist); % normalization
    if ii==1
        HhistMin = Hhist;
        HhistMax = Hhist;
    else
        HhistMin = min(HhistMin, Hhist);
        HhistMax = max(HhistMax, Hhist);
    end
end
close (h)