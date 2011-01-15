function [alpha,beta]=cnonlin2(evalstr,x,y,sy,p,delta,lambda);

% CNONLIN2 is a subsidiary function of CNONLIN. It must be called through CNONLIN.

m=length(p);
beta=-0.5*firstd(evalstr,x,y,sy,p,delta);
d2XiSq_da2=zeros(m);
for j=1:m;
  p(j)=p(j)+0.5*delta(j);
  dXiSqPlus=firstd(evalstr,x,y,sy,p,delta/1000);
  p(j)=p(j)-delta(j);
  dXiSqMinus=firstd(evalstr,x,y,sy,p,delta/1000);
  p(j)=p(j)+0.5*delta(j);
  d2XiSq_da2(:,j)=(dXiSqPlus-dXiSqMinus)/delta(j);
 end;
alpha=0.5*d2XiSq_da2;
for i=1:m
 alpha(i,i)=alpha(i,i)*(1+lambda);
end;
