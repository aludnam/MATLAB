% This computes FREM for different number of states the sources can get
% into (uniformly distributed over these states). 
int1_multi{1}=1;
int2_multi{1}=1;
int1_multi{2}=[0 1];
int2_multi{2}=[0 1];
int1_multi{3}=[0 .5 1];
int2_multi{3}=[0 .5 1];
int1_multi{4}=[0:.1:1];
int2_multi{4}=[0:.1:1];
clear('y_multi')

for mm=1:4
    
int1_vec= int1_multi{mm};
int2_vec= int2_multi{mm};
gFREMdemo
y_multi(:,mm)=y';
end

figure; 
plotdistanceDistributions

