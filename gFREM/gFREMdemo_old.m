showim=1;
x=[-15:.01:25];

% positions
l1=0;
l2=[0:.1:4, 7];

% sigma
sig1=1;
sig2=1;

% intensities
int1= 1;
int2_vec=[1];



for jj=1:length(int2_vec)
    int2=int2_vec(jj);
    for ii=1:length(l2)
        f1=makeGauss(x,l1,sig1);
        f2=makeGauss(x,l2(ii),sig2);
        y(ii,jj)=gFREMfunction(x,f1,f2, int1, int2,showim);
        if showim
            ylim([0,int1/(sqrt(2*pi)*sig1)+int2/(sqrt(2*pi)*sig2)])
            vline2(sig1*.61*2*pi/sqrt(2),'k--',{'Raleigh'}); % [Zhang et al., 2007]
            vline2(sig1*.47*2*pi/sqrt(2),'k--'); % [Zhang et al., 2007]
        end
        
    end
    y1(jj)=gFREMfunction(x,f1,f2,int1,0);
    y2(jj)=gFREMfunction(x,f1,f2,0,int2);
    l{jj}=num2str(['int [' num2str([int1 int2]) ']']);
end

plotdistance
