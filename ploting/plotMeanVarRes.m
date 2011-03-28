function [ms, ss, vs]=plotMeanVarRes(S)
% [ms, ss, vs]=plotMeanVarRes(S)
open ../S245/MeanVarQD_aplha17.fig
hold on

ms=mean(S,2)/10^3;
ss=std(S,[],2)/10^3;
vs=var(S,[],2)/10^6;
scatter(ms,vs,10^2,'dk','filled')
% for ii=1:size(S,1)
%     scatter(mean(S(ii,:))/1000, std(S(ii,:))/1000,1000,'xk')
% end