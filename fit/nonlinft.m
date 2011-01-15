function [pbest,perror,nchi2]=nonlinft(mfunc,x,y,sy,pt,v);
% function [pbest,perror,nchi2]=nonlinft(mfunc,x,y,sy,pt,v);
% 
% Levenberg-Marquardt non-linear regression.
% mfunc=name of the function file that calculates the function to be 
%    fitted to the set of data. This function cannot be complex; to fit
%    a complex function use cnonlin. mfunc is a string variable, so 
%    the name of the file must be put inside quotes, eg. 'model'.
%    The function file written to calculate your function should
%    be of the form:
%
%    function yfit=model(x,p)
%    yfit=p(1)+p(2)*x;
%    yfit=yfit(:);
%
%    The example above implements y=a+bx.
%    The vector x contains the values of the independent variable.
%    The vector p contains whatever parameters you are using. 
% 
% x=vector of the independent variable (e.g. Time).
% y=vector of the dependent variable (e.g. Counts recorded).
% sy=vector of standard error in y values.
% pt=initial estimate of parameters to be fitted. All parameters
%   must be non-zero, and are assumed to be real.
% v=vector indicating which parameters are to be varied (fitted), and 
%   which are to be held fixed. In most cases we will want all parameters 
%   to be varied, but in some circumstances it is useful to hold certain 
%   parameters fixed while the others are varied. v is a vector of the same
%   length as p, and should contain only ones and zeros; a one indicates that
%   the corresponding parameter should be varied (fitted), and a zero 
%   indicates that the corresponding parameter should be held constant.
%   eg. [1,0,1,1] would be used to keep the second parameter fixed while 
%   varying the rest (in the case of four parameters).
% 
% OUTPUT VARIABLES
% pbest=fitted parameters
% perror=standard error in fitted parameters
% nchi2=normalised chi-squared parameter (=chi-Squared/degrees of freedom). This
%    is expected to be close to 1 for a good fit.

% Written 4/12/95 by Michael Fleming, University Of Auckland

chiCut=0.01;  lambda=0.001;  stepSize=0.001;  chiOld=Inf;  
y=y(:);  x=x(:);  sy=sy(:);  pt=pt(:);  v=v(:);  maxIter=30;
m=length(pt); n=length(x);
 
evalstr=[mfunc '(x,['];      % Create string containing function call
count=0;
for i=1:m             
 if v(i)==1         % Extract vector p of parameters which will be varied
  count=count+1;    % from vector pt which contains all input parameters
  p(count)=pt(i);
  evalstr=[evalstr 'p(' int2str(count) ') '];
 else 
  evalstr=[evalstr num2str(pt(i)) ' '];
 end;
 if (pt(i)==0)   % Check all parameters are non-zero
  error('You have entered a zero parameter.All parameters must be nonzero');
 end;
end;
evalstr=[evalstr '])'];

dof=n-count;    p=p(:);   delta=p*stepSize;
% delta is the size of the small increment used for calculating numerical derivatives


% dof is the number of degrees of freedom


chiSqr=sum((abs(y-eval(evalstr))./sy).^2);  % Calc ChiSqr from initial parameters
if chiSqr/dof>5000;
 disp('You have made a bad choice of initial parameters');
end;

iter=0;
while (abs(chiOld-chiSqr)>chiCut)&(iter<maxIter);
 iter=iter+1;
 chiOld=chiSqr;
 [alpha,beta]=nonlin2(evalstr,x,y,sy,p,delta,lambda);
 if det(alpha)==0;
  error('No convergence - try a different set of parameters');
 end;
 dp=alpha\beta;    % Evaluate parameter increments
 p=p+dp;
 chiSqr=sum((abs(y-eval(evalstr))./sy).^2);
 while (chiSqr>(chiOld+chiCut));
  p=p-dp;
  iter=iter+1;  lambda=lambda*10;
  [alpha,beta]=nonlin2(evalstr,x,y,sy,p,delta,lambda);
  if det(alpha)==0;
   error('No convergence - try a different set of parameters');
  end;
  dp=alpha\beta;    % Evaluate parameter increments
  p=p+dp;
  chiSqr=sum((abs(y-eval(evalstr))./sy).^2);
 end;
 lambda=0.1*lambda;
end;

if iter==maxIter,
 disp('Maximum number of iterations exceeded - convergence not achieved');
end;

%xfine=linspace(x(1),x(n),201);    % Print graph of function
%yfitted=model(xfine,p);
%errorbar(x,y,sy);
% Calculate correlation coefficicent R^2
%f=eval(evalstr);   
%r=corrcoef(y,f);
%r2=r(1,2).^2;

% Calculate standard errors in parameters
[alpha,beta]=nonlin2(evalstr,x,y,sy,p,delta,0);	% Calculate final alpha matrix
sp=sqrt(diag(inv(alpha)));      % Evaluate standard errors in parameters.
disp('Parameters');
count=0;
for i=1:m
 if v(i)==1 
  count=count+1;
  pbest(i)=p(count); perror(i)=sp(count);
  disp(['p(', int2str(i), ') = ' sprintf('%5.3e',p(count)), '+-' sprintf('%5.3e',sp(count))]);
 else
  disp(['p(', int2str(i), ') = ' sprintf('%5.3e',pt(i)), '+-0']);
  pbest(i)=pt(i); perror(i)=0;
 end;
end;

format compact
%disp('Correlation Coefficient R^2');
%disp(r2);
disp('Normalised chi-squared value');
nchi2=chiSqr/dof;
disp(nchi2);
format;
