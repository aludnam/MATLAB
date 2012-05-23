function a=makeBars(numbers,sx,sy,scale,showim)
% makeBars creates images (sx X sy) of horizontal bars filled by ammount corresponding to numbers.
%
% a=makeBars(numbers,sx,sy,scale,showim)
% numbers - numbers to be represented
% sx,sy, dimensions of hte bar
% scale - if set to 1 (default) scales numbers between 0 and 1, else
% divides the numbers by their maximum
% showim - if set to 1 (default) displays images of bar into one line

if ~exist('scale','var')
    scale = 1;
end
if ~exist('showim','var')
    showim=1;
end

sn=length(numbers);
minn=min(numbers);
maxn=max(numbers);
if scale
    numsc=(numbers-minn)/(maxn-minn);
else 
    numsc=numbers/maxn; 
end
a=zeros(sy,sx,sn)+eps;
for ii=1:sn
    a(:,1:round(sx*numsc(ii)),ii)=1;
    a(:,1,ii)=1;
    a(:,end,ii)=0;
    
end

if showim
    figure; 
    imstiled(a,[],'gray',[],[1,sn])
end