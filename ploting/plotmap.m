function plotmap(X1,X2,plotwhat,plotname,zvec,limvec1,limvec2,scalecont,textONOFF)
% function
% plotmap(X1,X2,plotwhat,plotname,zvec,limvec1,limvec2,scalecont,textONOFF)
% zvec = [1 2 3 4];
% limvec1=[1:15];
% limvec2=[1:15];
% scalecont=20;
% plotname=xcovpeakr;
ihstr = [10 100 500 10000];

for ii=zvec
    figure
    contour(X1(limvec1,limvec2),X2(limvec1,limvec2),plotwhat(limvec1,limvec2,ii),scalecont, 'ShowText',textONOFF); 
    plottestNMF2
    
    SaveImageFULL([plotname '_iter' num2str(ihstr(ii))],'epf');
end