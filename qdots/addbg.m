function [w,h]=addbg(win, hin, bg)
% [w,h]=addbg(win, hin, bg)
% add the background component

sizew=size(win);
sizeh=size(hin);
w_bg = ones(sizew(1), 1)/sizew(1); %background component
w = [win, w_bg];

h_bg = bg*size(w,1)*ones(1, sizeh(2)); %background component
h=[hin; h_bg];
