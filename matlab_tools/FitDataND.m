% FitDataND : Fits a function to N-dimensional image data
% Synopsis: [fitted,params,myfunct] = FitData(AFunctString,InitParm,Data,maxiter,placement)
% AFunctString should use the variable 'x' as a vector and x{1}, x{2} to access its components, 
% and the values to fit: c(1),c(2),...
%
% InitParm : Vector of initial parameters c(..)
% Data : Experimental multidimensional (e.g. image) data to be fitted
% maxiter : maximum number of iterations in fit (default = 300)
% method : options are mean square error 'mse' and a fast version of i-divergence 'idiv', default 'mse'
% placement : position of zero coordinate in position vectors x{i}. 
%             Type 'help ramp' for explanation (default ='right')
%
%  Author : Rainer Heintzmann
%
% Example:
% a=noise(10*exp(-((xx(20,20)-2).^2+(yy(20,20)-1.2).^2)/20),'poisson') 
% [fitted,params,idiv,myfunct] = FitDataND('c(1)*exp(-((x{1}-c(2)).^2+(x{2}-c(3)).^2)/c(4))',[8 0 0 15],a,700,'idiv');
% fprintf('Intensity : %g, Position : (%g,%g), sigma^2: %g\n',params);
% fitted-a
% idiv

function [fitted,params,msevalue,myfunct] = FitDataND(AFunctString,InitParm,Data,maxiter,method,placement)

if nargin < 4
    maxiter = 700;
end

if nargin < 5
    method = 'mse';
end

if nargin < 6
    placement = 'right';
end

if strcmp(method,'mse')
    minfunc=@fiterror;
elseif strcmp(method,'idiv')
    minfunc=@fiterrorIdiv;
else
    tmp=['Unknown method: ' method '. Use ''mse'' or ''idiv'''];
    error(tmp);
end

myfunct = inline(AFunctString,'c','x');
options=optimset('Display','notify','TolX',10^-4,'TolFun',10^-8,'MaxIter',maxiter);

for d=1:size(size(Data),2)
    Coords{d}=ramp(size(Data),d,placement);
end

[params,msevalue]=fminsearch(minfunc,InitParm,options,myfunct,Data,1,Coords);

fitted = myfunct(params,Coords);


% here is a function defining the error estimate
function msevalue=fiterror(params,myfunct,Data,wght,Coords)
  for d=1:size(size(Data),2)
    simdata=myfunct(params,Coords);
  end
  msevalue = mean(wght.* ((Data-simdata)).^2)/sum(wght);
%msevalue

function fidivval=fiterrorIdiv(params,myfunct,Data,wght,Coords)  % This is the fast i-divergence function
  for d=1:size(size(Data),2)
    simdata=myfunct(params,Coords);
  end
  fidivval=mean(simdata-Data .* log(simdata));  % fast version omitting constants
