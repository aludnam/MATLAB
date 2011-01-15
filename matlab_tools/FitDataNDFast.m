% FitDataNDFast : Fits several ND Gaussians to N-dimensional image data
% Synopsis: [fitted,params,myfunct] = FitDataNDFast(InitParm,Data,numdims,maxiter,method)
% AFunctString should use the variable 'x' as a vector and x{1}, x{2} to access its components, 
% and the values to fit: c(1),c(2),...
%
% InitParm : Vector of initial parameters first row is global parameters,
% meaning is as follows: [global background, global width, intensity, posy, posx, intensity2, posy2, posx2...]
% other rows are local parameters
% Data : Experimental multidimensional (e.g. image) data to be fitted
% maxiter : maximum number of iterations in fit (default = 300)
% method : figure of merit to use for fitting 'mes', 'idiv' or 'fidiv'
%
% Example:
% a=noise(51.2*exp(-((xx(20,20)-2).^2+(yy(20,20)-1.2).^2)/20)+33,'poisson')
% [params,res,fitted]=FitDataNDFast([30 15 40 2 2],a,2,300,'idiv')
% fitted is the residual image (depending on method)
% params contains the result of the fit

function [params,res,fitted,residual] = FitDataNDFast(InitParm,Data,numdims,maxiter,mymethod)

if nargin < 3
    numdims = 2;
end
if nargin < 4
    maxiter = 3000;
end

if nargin < 5
    mymethod = 'mse';
end

%if nargin < 5
%    placement = 'right';
%end

%fixedparams=[0 20];
MultiGaussMSE(double(Data),mymethod,numdims)

[res, fitted] = MultiGaussMSE(InitParm');
%params=fixedparams;


if (0)  % old method
options=optimset('Display','notify','TolX',10^-4,'TolFun',10^-8,'MaxIter',maxiter);
[params,msevalue]=fminsearch(@MultiGaussMSE,InitParm',options);
else    % new method with smarter amd faster optimisation routine from http://www.cs.ubc.ca/~schmidtm
options=struct('Display','off','notify',1,'numDiff',1,'TolX',10^-4,'TolFun',10^-8,'MaxIter',maxiter);
[params,msevalue,moreinfo]=minFunc(@MultiGaussMSE,InitParm',options);
end

[res, fitted,residual] = MultiGaussMSE(params);
params = params';
fitted=dip_image(fitted);
residual=dip_image(residual);
