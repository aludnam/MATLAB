function a=makeBars(numbers,sx,sy,scale,showim,maxnumber)
% makeBars creates images (sx X sy) of horizontal bars filled by ammount corresponding to numbers.
%
% a=makeBars(numbers,sx,sy,scale,showim,maxnumber)
% numbers - numbers to be represented
% sx,sy, dimensions of hte bar
% scale - if set to 1 (default) scales numbers between 0 and 1 (unless there is maxnumber), else
% divides the numbers by their maximum
% showim - if set to 1 (default) displays images of bar into one line
% maxnumber - maximum number to which the rest is scaled

if ~exist('scale','var')
    scale = 1;
end
if ~exist('showim','var')
    showim=1;
end
if ~exist('maxnumber','var')
    maxnumber=0;
end
if maxnumber
    scale=1;
end

sn=length(numbers);
minn=min(numbers);
maxn=max(numbers);
if scale
    if ~maxnumber
        numsc=(numbers-minn)/(maxn-minn);
    else
        numsc=numbers/maxnumber;
    end
else 
    numsc=numbers/maxn; 
end
if sum(numsc>1)
    warning('Problem with scaling the bars!')
    numsc=min(numsc,1); % when the number is over one it is clipped to 1.
end
    
a=zeros(sy,sx,sn)+eps;
for ii=1:sn
    a(:,1:ceil(sx*numsc(ii)),ii)=1; 
    a(:,1,ii)=1;
    a(:,end,ii)=0;
    
end

if showim
    figure; 
    imstiled(a,[],'gray',[],[1,sn])
end