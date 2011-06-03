function [pint, int_out]=generateDistribution(int_vec,tau, probfunction, correctIntensity)
% [pint, int_out]=generateDistribution(int_vec,tau, probfunction, correctIntensity)
% Generates distributions over states for gFREM.
if ~exist('correctIntensity','var')
    correctIntensity = 0; 
end
int1_vec=int_vec(1,:);
int2_vec=int_vec(2,:);
tau1=tau(1,:);
tau2=tau(2,:);
lint1=length(int1_vec);
lint2=length(int2_vec);
meanint1 = .5; %mean intensity of the source 1
meanint2 = .5;
if size(int_vec,2) > 1
    switch probfunction
        case 'uniform'
            pint1=1/lint1*ones(1,lint1);
            pint2=1/lint2*ones(1,lint2);
        case 'exponential'
            aa1=tau1;
            aa2=tau2;
            pint1=1/sum(exp(-int1_vec/aa1))*exp(-int1_vec/aa1);
            pint2=1/sum(exp(-int2_vec/aa2))*exp(-int2_vec/aa2);
            %         figure(20); hold on; plot(int1_vec,pint1,'o-');
        case 'binomial'
            pint1=[tau1 1-tau1];
            pint2=[tau2 1-tau2];
            
    end
else 
    pint1=1;
    pint2=1;
end

if correctIntensity % correct intenstity so that the mean is kept constant. 
    fprintf('Intensity vector adusted by a factor %g to keep the mean constant.\n',1/(int1_vec*pint1'))
    int1_vec=meanint1*int1_vec/(int1_vec*pint1'); %to keep the mean constant
    int2_vec=meanint2*int2_vec/(int2_vec*pint2');
end
pint=cat(1,pint1,pint2);
int_out=cat(1,int1_vec,int2_vec);
