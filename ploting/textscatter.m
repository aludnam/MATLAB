function textscatter(varargin)
% textscatter(cx,cy,markersize,markerstyle,textcolor, textoffsetx, textoffsety)
% Uses scatter plot and each marker is labeled by the number of the row in
% the [cx, cy] matrix. Text is offste by [textoffsetx, textoffsety] from
% the markers. 
% Examples: 
% textscatter(100*rand(1,30), 100*rand(1,30))
% textscatter(100*rand(1,30), 100*rand(1,30),[],'d')
% textscatter(100*rand(1,30), 100*rand(1,30),100,'gx','r', 2,3)

args = varargin;    
if nargin<5 % textcolor 
%     args{5} = [0 0 0]; % black   
    args{5} = [1 0 0]; % red   
end
if nargin<6 % textoffset
    args{6}= .5;
    args{7}= .5;
end

cx = args{1}; 
cy = args{2};
textcolor = args{5};
textoffsetx = args{6};
textoffsety = args{7};

scatter(args{1:min(nargin,4)})
text(cx+textoffsetx, cy+textoffsety, num2str([1:length(cx)]' ),'color' , textcolor)