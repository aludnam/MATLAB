function K=fitTwoLines(x,y, showimage) 
% K=fitTwoLines(x,y, showimage)
% This fits two adjecent lines into data [x,y] which are partitioed as
% x=[x1,x2] and y=[y1,y2]. It determines K which gives best fit for
% partition x1=x(1:K), x2=x(K:end) ...
% Used for estimation of the 'kink' in the graph.  
if ~exist('showimage','var')
    showimage = 0; 
end
n=length(y);

for ind_includeMidPoint=0:1 % when 0 then middle point is included (shared by both lines)
    index =0;
    for ii=1:n-ind_includeMidPoint
        index = index + 1;
        x1=x(1:ii);
        y1=y(1:ii);
        x2=x(ii+ind_includeMidPoint:n);
        y2=y(ii+ind_includeMidPoint:n);
        
        [pbest,perror,nchi]=nonlinft('linearfunction' ,x1,y1,ones(size(x1)),[-1 1e-5],[1 1]);
        pbest1(index,:,ind_includeMidPoint+1)=pbest;
        perror1(index,:,ind_includeMidPoint+1)=perror;
        nchi1(index,ind_includeMidPoint+1)=nchi;
        [pbest,perror,nchi]=nonlinft('linearfunction' ,x2,y2,ones(size(x2)),[-1 1e-5],[1 1]);
        pbest2(index,:,ind_includeMidPoint+1)=pbest;
        perror2(index,:,ind_includeMidPoint+1)=perror;
        nchi2(index,ind_includeMidPoint+1)=nchi;
    end
end
nn=squeeze((sqrt(perror2(:,1,:).^2+perror2(:,2,:).^2)));
objectiveFunction = (nchi1.^2+nchi2.^2);
mo=min(objectiveFunction,[],2); % minimum the two either including (1) or not (2) the middle point 
[del, ixmo]=min(objectiveFunction./(1./nn),[],2); %whether to choose (1) or (2). This is weighted by error of the estimation...
% [m,mx]=min(objectiveFunction);
% [m,mx]=min(mo)+ixmo-1;
% K=x(mx+1);
[m,mx]=min(mo(1:end-1));
K=x(mx+ixmo(mx)-1);

if showimage
    figure;
    plot(x(1:end), mo)
    vline2(K,'r','K');    
    xlabel ('x')
    ylabel ('objective function')
    
    figure;
    plot(x,y,'-o')
    hold on
    plot(x(1:mx), pbest1(mx,1,ixmo(mx))*x(1:mx)+pbest1(mx,2,ixmo(mx)),'g')
    plot(x(mx+ixmo(mx)-1:n), pbest2(mx,1, ixmo(mx))*x(mx+ixmo(mx)-1:n)+pbest2(mx,2,ixmo(mx)),'r')
    xlabel('x')
    ylabel('y')
end