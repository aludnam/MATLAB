% FitNlinData : Fits a function to the data
% Synopsis: [fitted,params,myfunct,msevalue,redchi2,limit] = FitNlinData(AFunctString,InitParm,DataY,DataX,maxiter)
% DataX: default is 1:size(DataY)
% maxiter: default 700 
% method : options are mean square error 'mse' and a fast version of i-divergence 'idiv', default 'mse'
% AFunctString should use the variable 'x', and the values to fit: c(1),c(2),...
% msevalue : mean square error
% redchi2 : reduced chi^2
% limit : approximation to the 5% limit of redchi2 to estimate whether to
%         reject this function as a good fit.
% Autor : Rainer Heintzmann
%
% Example:
% x= 1:100;
% y=noise(2.5*x.^2-212*x+5000,'gaussian',300);
% [fitted,params,myfunct,mse] = FitNlinData('c(1)*x.^2+c(2)*x+c(3)',[1 1 1],y);
% plot(x,y)
% hold on
% plot(x,fitted,'red')
% fprintf(' Function is : %g*x^2+%g*x+%g\n',params)

function [fitted,params,myfunct,msevalue,redchi2,limit] = FitNlinData(AFunctString,InitParm,DataY,DataX,maxiter,method)
if (size(DataY,1) > size(DataY,2))
    DataY=DataY';
end

if (nargin < 4)
    DataX=1:size(DataY);
end

if (size(DataX,1) > size(DataX,2))
    DataX=DataX';
end

if nargin < 5
    maxiter = 700;
end

if nargin < 6
    method = 'mse';
end


myfunct = inline(AFunctString,'c','x');
if strcmp(method,'mse')
    minfunc=@fiterror;
elseif strcmp(method,'idiv')
    minfunc=@fiterrorIdiv;
else
    tmp=['Unknown method: ' method '. Use ''mse'' or ''idiv'''];
    error(tmp);
end

options=optimset('Display','notify','TolX',10^-4,'TolFun',10^-8,'MaxIter',maxiter);

[params,msevalue]=fminsearch(minfunc,InitParm,options,myfunct,DataX,DataY,1);

%params = nlinfit(DataX,DataY,myfunct,InitParm);
fitted = myfunct(params,DataX);

resid = DataY-fitted;
dresid = resid(1:end-1)-resid(2:end);
stddev = sqrt(var(dresid)/2);
redchi2 = sum((DataY-fitted).^2 ./ (stddev^2))/(size(fitted,2) - size(InitParm,2));
limit = 1 + (2*sqrt(2))/sqrt(size(fitted,2)- size(InitParm,2));

% here is a function defining the error estimate
function msevalue=fiterror(params,myfunct,DataX,DataY,wght)
  simdata=myfunct(params,DataX);
  msevalue= mean(wght.* ((DataY-simdata)).^2)/sum(wght);
  
function fidivval=fiterrorIdiv(params,myfunct,DataX,DataY,wght)  % This is the fast i-divergence function
  simdata=myfunct(params,DataX);
  fidivval=mean(simdata-DataY .* log(simdata));  % fast version omitting constants
