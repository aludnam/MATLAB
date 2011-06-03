showim=0;
x=[-15:.01:25];

% positions
l1=0;
l2=[0:.1:4,5.5, 7];

% sigma
sig1=1;
sig2=1;

% intensities
% int1_vec= [0:.1:1];
% int2_vec=[0:.1:1];
lint1=length(int1_vec);
lint2=length(int2_vec);
% distribution of the intensities

probfunction = 'exponential';
switch probfunction
    case 'uniform'
        pint1=1/lint1*ones(1,lint1);
        pint2=1/lint2*ones(1,lint2);
    case 'exponential'        
        meanint1 = .5; %mean intensity of the source 1
        meanint2 = .5;
        aa1=tau1;
        aa2=tau2;
        pint1=1/sum(exp(-int1_vec/aa1))*exp(-int1_vec/aa1);
        pint2=1/sum(exp(-int2_vec/aa2))*exp(-int2_vec/aa2);
        fprintf('Intensity vector adusted by a factor %g to keep the mean constant.\n',1/(int1_vec*pint1'))
        int1_vec=meanint1*int1_vec/(int1_vec*pint1'); %to keep the mean constant
        int2_vec=meanint2*int2_vec/(int2_vec*pint2'); 
        figure(20); hold on; plot(int1_vec,pint1,'o-');
    case 'gamma'
        
end


y=zeros(1,length(l2));
for ind_dist=1:length(l2)
    for ind_int1=1:lint1
        int1=int1_vec(ind_int1);
        for ind_int2=1:lint2
            int2=int2_vec(ind_int2);            
            f1=makeGauss(x,l1,sig1);
            f2=makeGauss(x,l2(ind_dist),sig2);
            %             y(ind_dist,ind_int2)=gFREMfunction(x,f1,f2, int1, int2,showim);            
            if and(int1==0, int2==0)
                y(ind_dist)=0;
            else
                y(ind_dist)=y(ind_dist)+pint1(ind_int1)*pint2(ind_int2)*gFREMfunction(x,f1,f2, int1, int2,showim);
            end
            if showim
                ylim([0,int1/(sqrt(2*pi)*sig1)+int2/(sqrt(2*pi)*sig2)])
                vline2(sig1*.61*2*pi/sqrt(2),'k--',{'Raleigh'}); % [Zhang et al., 2007]
                vline2(sig1*.47*2*pi/sqrt(2),'k--'); % [Zhang et al., 2007]
            end
            
        end
        y1(ind_int2)=gFREMfunction(x,f1,f2,int1,0);
        y2(ind_int2)=gFREMfunction(x,f1,f2,0,int2);
%         l{jj}=num2str(['int [' num2str([int1 int2]) ']']);
    end
end
l=[];
% plotdistance
