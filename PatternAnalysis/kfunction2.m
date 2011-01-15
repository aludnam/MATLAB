function K = kfunction2(data1, data2 , xK,box,method)
% KFUNCTION2 calculates Ripleys K function 
% K = kfunction2(data1, data2, xK, box, method) - returns vector K containing value
% of Ripley's K-function of data1 with respect to data2 in the distances in xK.
% data1 - N1-by-2 vector where N is number of datapoints. Each row
% corresponds to x and y coordinates of each datapoint
% data2 - N2-by-2 vector where N is number of datapoints. Each row
% corresponds to x and y coordinates of each datapoint
% xK - corresponds to the distances where K function should be computed.
% K is the same size as xK...
% box - rectangular boudnary of the data: box = [xlim1, xlim2, ylim1,
% ylim2]
% method - switch between edge correction. If method=0, no edge correction
% is applied. If method=1, datapoint is used for estimation of K(h) only if
% it is at least h units away from the box


if nargin<4 method=1; end
[N1,k1] = size(data1);
[N2,k2] = size(data2);
if or(k1~=2, k2~=2) error('data1 must have two columns'); end

%%%%%%%%!!!!!!!!!!! upravit
rbox1 = min([   data1(:,1)'-box(1);
                box(2)-data1(:,1)';
                data1(:,2)'-box(3);
                box(4)-data1(:,2)']);

% rbox1 is the nearest distance of each datapoint data1 to the box

data = [data1; data2];
DISTA = squareform(pdist(data,'euclidean'));

DIST = DISTA(N1+1:end, 1:N1); %distances between data1 and data2

% DIST = sort(DIST); %sorts along columns only

if method==1 % edge correction...
K = zeros(length(xK),1);
Nk = length(K);
wb = waitbar(0,'Computing Ripley''s K-function...');
for k=1:Nk
    waitbar(k/Nk,wb);    
    I = find(rbox1>=xK(k));
    if ~isempty(I)
        K(k) = sum(sum(DIST(1:end,I)<=xK(k)))/length(I);
    end
end
close (wb);

elseif method==0 % no correction
    K = zeros(length(xK),1);
    for k=1:length(K)
        K(k) = sum(sum(DIST(1:end,:)<=xK(k)))/N1;
    end
end

% fprintf ('\b'); fprintf ('\b'); fprintf ('\b'); fprintf ('\b');
% fprintf ('100%%\n');
A = (box(2)-box(1))*(box(4)-box(3));
lambda2 = N2/A;
K = K/(lambda2);