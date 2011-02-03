function [wout, hout] = nomalize_wh(win,hin)
% [wout, hout] = nomalize_wh(win,hin,k)

sumw = sum(win,1);
wout = win./repmat(sumw,size(win,1),1);                %normalization of each component
hout = hin.*repmat(sumw',1,size(hin,2));    %to keep the multiplication equal
