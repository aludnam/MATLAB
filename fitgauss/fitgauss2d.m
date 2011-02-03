function [x_mu, y_mu, sig, differ] = fitgauss2d(M,sigfix, showfig)
% [x_mu, y_mu, sig, differ] = fitgauss2d(M,sigfix, showfig)
% x_mu, y_mu : coordiantes of the fitted mean
% sig : fitted std
% differ : D -devergence between data and figure generated form gaussian


M = abs(M)+1e-9;
mM = max(max(M));
[xm, ym] = find(mM==M);
M = M/mM;
if ~exist('sigfix','var')
    sigfix = [];
end
if isempty(sigfix)
    sguess = sum(sum(M))/(size(M,1)*size(M,2))*size(M,1);
    x = fminsearch(@(x) difference(x,M),[ym(1),xm(1),sguess]);
    x_mu = x(1);
    y_mu = x(2);
    sig = x(3);
else
    x = fminsearch(@(x) difference_sigfix(x,M,sigfix),[ym,xm]);
    x_mu = x(1);
    y_mu = x(2);
    sig = sigfix;
end

differ = difference(x,M);

if exist('showfig','var')
    if showfig == 1
        
%         ims(gauss2d(size(M), [x_mu y_mu], sig,1));
        ims(M, 'gray');
        hold on
        scatter (x_mu, y_mu, 80 , 'xr');
        [x,y,z] = cylinder(sig,200);
        plot(x(1,:) + x_mu, y(1,:) + y_mu,'r')
    end
end