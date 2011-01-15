function [pbest, perror, nchi2, options] = fminerr(fitfcn,x,options,P1,P2,P3)
%fminerr   Fits non-linear curves to experimental data.
%	[pbest, perror,nchi2]=fminerr('F',p,OPTIONS,xdata,ydata,sy) attempts to return a 
%  	vector pbest which is the best fit set of parameters for fitting the function 
%       'F' to a set of data near the starting guess p. The standard errors in the 
%       parameters are returned in the vector perror. 'F' is a string containing the
%	name of the objective function to be fitted.  F(x) should be a
%	scalar valued function of a vector variable. The vectors xdata and ydata
%       contain the experimental data, and sy is a vector of standard errors in
%       ydata.
%	If OPTIONS(1) is nonzero, intermediate steps in the solution are
%	displayed; the default is OPTIONS(1) = 0.  OPTIONS(2) is the termination
%	tolerance for x; the default is 1.e-4.  OPTIONS(3) is the termination
%	tolerance for F(x); the default is 1.e-4.  OPTIONS(14) is the maximum
%	number of steps; the default is OPTIONS(14) = 500.  The other components
%	of OPTIONS are not used as input control parameters by FMINERR.  For more
%	information, see FOPTIONS.
%
%       OUTPUT: 
%       pbest=best fit parameters
%       perror=standard error in pbest
%       nchi2=normalised chi-squared value
%
%	FMINERR uses a Simplex search method.
%
%	See also FMINS. 

%	Reference: J. E. Dennis, Jr. and D. J. Woods, New Computing
%	Environments: Microcomputers in Large-Scale Computing,
%	edited by A. Wouk, SIAM, 1987, pp. 116-122.

%	C. Moler, 8-19-86
%	Revised Andy Grace, 6-22-90, 1-17-92 CBM,  10-5-93 AFP
%	Copyright (c) 1984-94 by The MathWorks, Inc.

%       Revised by Mike Fleming 29/11/95 for non-linear data fitting with
%       parameter error estimations.

if nargin<3, options = []; end
options = foptions(options);
prnt = options(1);
tol = options(2);
tol2 = options(3);

evalstr = [fitfcn,'(P1,x)'];

n = prod(size(x));
if (~options(14)) 
    options(14) = 200*n; 
end

P1=P1(:); P2=P2(:); P3=P3(:);  % Ensure column vectors

% Set up a simplex near the initial guess.
xin = x(:); % Force xin to be a column vector
v = xin; 
x(:) = v; fv = sum((abs(P2-eval(evalstr))./P3).^2); 

% Following improvement suggested by L.Pfeffer at Stanford
usual_delta = 0.05;             % 5 percent deltas for non-zero terms
zero_term_delta = 0.00025;      % Even smaller delta for zero elements of x
for j = 1:n
   y = xin;
   if y(j) ~= 0
      y(j) = (1 + usual_delta)*y(j);
   else
      y(j) = zero_term_delta;
   end
   v = [v y];
   x(:) = y; f = sum((abs(P2-eval(evalstr))./P3).^2);
   fv = [fv  f];
end
[fv,j] = sort(fv);
v = v(:,j);


cnt = n+1;
if prnt
   clc
   format compact
   format short e
   home
   cnt
   disp('initial ')
   disp(' ')
   v
   f
end

alpha = 1;  beta = 1/2;  gamma = 2;
[n,np1] = size(v);
onesn = ones(1,n); 
ot = 2:n+1;
on = 1:n;

% Iterate until the diameter of the simplex is less than tol.
while cnt < options(14)
    if max(max(abs(v(:,ot)-v(:,onesn)))) <= tol & ...
           max(abs(fv(1)-fv(ot))) <= tol2
        break
    end

    % One step of the Nelder-Mead simplex algorithm

    vbar = (sum(v(:,on)')/n)';
    vr = (1 + alpha)*vbar - alpha*v(:,n+1);
    x(:) = vr;
    fr = sum((abs(P2-eval(evalstr))./P3).^2); 
    cnt = cnt + 1; 
    vk = vr;  fk = fr; how = 'reflect ';
    if fr < fv(n)
        if fr < fv(1)
            ve = gamma*vr + (1-gamma)*vbar;
            x(:) = ve;
            fe = sum((abs(P2-eval(evalstr))./P3).^2);
            cnt = cnt + 1;
            if fe < fv(1)
                vk = ve; fk = fe;
                how = 'expand  ';
            end
        end
    else
        vt = v(:,n+1); ft = fv(n+1);
        if fr < ft
            vt = vr; ft = fr;
        end
        vc = beta*vt + (1-beta)*vbar;
        x(:) = vc;
        fc = sum((abs(P2-eval(evalstr))./P3).^2); 
        cnt = cnt + 1;
        if fc < fv(n)
            vk = vc; fk = fc;
            how = 'contract';
        else
            for j = 2:n
                v(:,j) = (v(:,1) + v(:,j))/2;
                x(:) = v(:,j);
                fv(j) = sum((abs(P2-eval(evalstr))./P3).^2); 
            end
        cnt = cnt + n-1;
        vk = (v(:,1) + v(:,n+1))/2;
        x(:) = vk;
        fk = sum((abs(P2-eval(evalstr))./P3).^2); 
        cnt = cnt + 1;
        how = 'shrink  ';
        end
    end
    v(:,n+1) = vk;
    fv(n+1) = fk;
    [fv,j] = sort(fv);
    v = v(:,j);

    if prnt
        home
        cnt
        disp(how)
        disp(' ')
        v
        fv
    end
end
x(:) = v(:,1);
if prnt, format, end
options(10)=cnt;
options(8)=min(fv); 
if cnt==options(14) 
    if options(1) >= 0
        disp(['Warning: Maximum number of iterations (', ...
               int2str(options(14)),') has been exceeded']);
        disp( '         (increase OPTIONS(14)).')
    end
end

P1=P1(:);
P2=P2(:);
P3=P3(:);
x=x(:);
m=length(x);
delta=0.001*x;

% Calculate error matrix using second derivative of chi-squared
d2XiSq_da2=zeros(m);
for j=1:m;
  x(j)=x(j)+0.5*delta(j);
  for i=1:m;
   x(i)=x(i)+0.5*delta(i);
   XiSqPlus=sum((abs(P2-eval(evalstr))./P3).^2);
   x(i)=x(i)-delta(i);
   XiSqMinus=sum((abs(P2-eval(evalstr))./P3).^2);
   x(i)=x(i)+0.5*delta(i);
   dXiSq_da(i)=(XiSqPlus-XiSqMinus)/delta(i);
  end;
  dXiSqPlus=dXiSq_da(:);
  x(j)=x(j)-delta(j);
  for i=1:m;
   x(i)=x(i)+0.5*delta(i);
   XiSqPlus=sum((abs(P2-eval(evalstr))./P3).^2);
   x(i)=x(i)-delta(i);
   XiSqMinus=sum((abs(P2-eval(evalstr))./P3).^2);
   x(i)=x(i)+0.5*delta(i);
   dXiSq_da(i)=(XiSqPlus-XiSqMinus)/delta(i);
  end;
  dXiSqMinus=dXiSq_da(:);
  x(j)=x(j)+0.5*delta(j);
  d2XiSq_da2(:,j)=(dXiSqPlus-dXiSqMinus)/delta(j);
 end;
alpha=0.5*d2XiSq_da2;
perror=sqrt(diag(inv(alpha)));
pbest=x;

disp('  Parameters   Standard Errors');
format short e
disp([pbest,perror]);
format
format compact
disp('Normalised chi-squared value');
dof=length(P1)-m;
nchi2=sum((abs(P2-eval(evalstr))./P3).^2)/dof;
disp(nchi2);
format;
