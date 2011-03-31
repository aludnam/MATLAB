function [x1in,x2in]=initx1x2(Hinit,Winit,x1,x2,x1rand,x2rand)
% [x1in,x2in]=initx1x2(Hinit,Winit,x1,x2,x1rand,x2rand)
% initialization of the W and H for fitting
if strcmp(Hinit,'trueH')
    x1in=log(x1);
elseif strcmp(Hinit,'wrongH')
    x1in=log(x1rand);
else
    error('wrong name')
end

if strcmp(Winit,'trueW')
    x2in = x2;
elseif strcmp(Winit,'wrongW')
    x2in = x2rand;
else
    error('wrong name')
end