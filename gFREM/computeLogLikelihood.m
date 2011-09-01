function [logl, c1, c2, p] = computeLogLikelihood(d,bg)
% [logl, c1, c2, p] = computeLogLikelihood(d,bg)
p.dimensionality =2; % Number of diensions (1 or 2) of the PSF.
p.d=d;
p.bg=bg;
p.lambda = 655; %nm
p.NA = 1.2;
p.pixelsize = 106; %nm
p.sig1=sqrt(2)/2/pi*p.lambda/p.NA/p.pixelsize; %[Zhang 2007]
p.sig2=p.sig1;
p.int1=750;
p.int2=750;

% x=-25:.1:25;
p.stepCoord = 0.1;
% x1=-7:stepCoord:16;
% boundnary = 25; 
boundnary = 10; 
x1=-boundnary:p.stepCoord:boundnary;
y1=-boundnary:p.stepCoord:boundnary;

if p.dimensionality ==1
    % This is for 1D evaluation:
    x = x1;    
    % x=-7:.5:16;
    % x=-10:.01:20;
elseif p.dimensionality ==2
    % This is for 2D evaluation:
    [xx,yy]=meshgrid(x1,y1);
    x=cat(3,xx,yy);
end

p.c1_true = p.d/2; %[0 .5 1 1.5 2 3 4];
p.c2_true = -p.c1_true;
q=8;
qstep = .5;
% c1=[p.c1_true-q:.2:p.c1_true+q];
% c2=[p.c2_true-q:.2:p.c2_true+q];
c1=[-q:qstep:q];
c2=[-q:qstep:q];
if ndims(x) == 2 %1D vector
    f1_true = makeGauss(x,p.c1_true,p.sig1);                    % creates PSF (gauss approx -> !!! Different to simulationtools/makegauss.m !!!)
    f2_true = makeGauss(x,p.c2_true,p.sig2);
else
    f1_2D=makeGauss2D(x,p.c1_true,p.sig1);
    f2_2D=makeGauss2D(x,p.c2_true,p.sig2);
    f1_true = reshape(f1_2D,1,numel(f1_2D));                    % making 1D vector by concatenating 2D array
    f2_true = reshape(f2_2D,1,numel(f2_2D));
end


l_true = p.int1*f1_true+p.int2*f2_true+p.bg;
logl=zeros(length(c1),length(c2));

for ii=1:length(c1)
    % creates PSF (gauss approx)
    if ndims(x) == 2 %1D vector
        f1=makeGauss(x,c1(ii),p.sig1);
    else
        f1_2D=makeGauss2D(x,c1(ii),p.sig1);
        f1 = reshape(f1_2D,1,numel(f1_2D));                    % making 1D vector by concatenating 2D array
    end
    
    for jj=1:length(c2)
        if ndims(x) == 2 %1D vector
            f2=makeGauss(x,c2(jj),p.sig2);            % creates PSF (gauss approx)
        else
            f2_2D=makeGauss2D(x,c2(jj),p.sig2);
            f2 = reshape(f2_2D,1,numel(f2_2D));                    % making 1D vector by concatenating 2D array
        end
        
        
        l=p.int1*f1+p.int2*f2+p.bg;
        logl(ii,jj)=sum(l_true.*log(l)-l-factorialapprox(l_true));
    end
end
