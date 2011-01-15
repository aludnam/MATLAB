function [xout, yout]=changerelativedistance(x,y,coef, offset_vec)
% [xout, yout]=changerelativedistance(x,y,coef, offset_vec)
% Stretches (or enlarges) distances between x and y from their mean by coef
% ratio and shifts them by offset_vec
% default offset_vec: 0;

if ~exist('offset_vec', 'var')
    offset_vec=[0 0];
end
mvx=mean(x);
mvy=mean(y);
mvxdistout=coef*(x-mvx);
mvydistout=coef*(y-mvy);
xout=mvx+mvxdistout+offset_vec(1);
yout=mvy+mvydistout+offset_vec(2);
