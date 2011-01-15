function [alpha,beta]=nonlin2(evalstr,x,y,sy,p,delta,lambda);

% nonlin2 is a subsidiary function of nonlinft. It must be called through nonlinft.

m=length(p);
n=length(x);
% Compute beta=-0.5*grad(chiSqr)
% Numerical partial derivatives of chi-squared are taken in an interval around point p
for i=1:m;
 p(i)=p(i)+0.5*delta(i);
 XiSqPlus=sum((abs(y-eval(evalstr))./sy).^2);
 p(i)=p(i)-delta(i);
 XiSqMinus=sum((abs(y-eval(evalstr))./sy).^2);
 p(i)=p(i)+0.5*delta(i);
 dXiSq_da(i)=(XiSqPlus-XiSqMinus)/delta(i);
end;
beta=-0.5*dXiSq_da;
beta=beta(:);

%make alpha
alpha=zeros(m);

y0=eval(evalstr);
for i=1:m;
 p(i)=p(i)+delta(i);
 y1=eval(evalstr);
 p(i)=p(i)-delta(i);
 dYda(:,i)=(y1-y0)./(sy*delta(i));
 end;
for j=1:n
 jdYda=dYda(j,:);
 alpha=alpha+jdYda'*jdYda;
end;

for i=1:m
 alpha(i,i)=alpha(i,i)*(1+lambda);
end;
