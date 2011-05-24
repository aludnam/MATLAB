jitter = .2;
ysp = 2; 
ymaxsp = 4+2;
xsp = 0.5;
xmaxsp = 5;
yfix = [2:ysp:ymaxsp]';
xfix = [0:xsp:xmaxsp]';
nxfix= length(xfix);
nyfix= length(yfix);

coord=zeros(nxfix*nyfix,2);
indstart = 1; 
indend= nxfix;
for ii=1:length(yfix)
    coord(indstart:indend,1)=xfix+jitter*rand(nxfix,1);
    coord(indstart:indend,2)=repmat(yfix(ii),nxfix,1)+jitter*rand(nxfix,1);
    indstart=indend+1;
    indend=indstart+nxfix-1;
end


jitter = .2;
ysp = .5; 
ymaxsp = 8;
xsp = 2;
xmaxsp = 3.5;
yfix = [0:ysp:ymaxsp]';
xfix = [1.5:xsp:xmaxsp]';
nxfix= length(xfix);
nyfix= length(yfix);

coord2=zeros(nxfix*nyfix,2);
indstart = 1; 
indend= nxfix;
for ii=1:length(yfix)
    coord2(indstart:indend,1)=xfix+jitter*rand(nxfix,1);
    coord2(indstart:indend,2)=repmat(yfix(ii),nxfix,1)+jitter*rand(nxfix,1);
    indstart=indend+1;
    indend=indstart+nxfix-1;
end