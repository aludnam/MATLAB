function plotMeanVarRes(S)
open ../S245/MeanVarQD_aplha17.fig
hold on

for ii=1:size(S,1)
    scatter(mean(S(ii,:))/1000, std(S(ii,:))/1000,1000,'xk')
end