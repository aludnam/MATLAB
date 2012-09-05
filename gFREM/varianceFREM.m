% This computes FREM for different number of states the sources can get
% into (uniformly distributed over these states).
dimensionality =2; % Number of dimensions (1 or 2) of the PSF.
savethis =0;
p.offset = 100     ;

int1_multi{1}=[2000];       % intensity of the source one
% int1_multi{2}=[2000];      % this allows comparison of sources with different intensity
% int1_multi{3}=[6000];
% int1_multi{4}=[6000];

int2_multi = int1_multi;    % intensity of the source two

% positions of sources
l1=0;
% l2=0:.2:10;
% l2=0:.05:1;
l2 = .4;

% coordinates of the images (pixelised version needs to do some binnign)
stepCoord = 0.2;
% x1=-7:stepCoord:16;
x1=-5:stepCoord:6;
y1=-4:stepCoord:4;

if dimensionality ==1
    % This is for 1D evaluation:
    x = x1;    
    % x=-7:.5:16;
    % x=-10:.01:20;
elseif dimensionality ==2
    % This is for 2D evaluation:
    [xx,yy]=meshgrid(x1,y1);
    x=cat(3,xx,yy);
end
fprintf('Computaiton for %gD PSF.\n',dimensionality);

p.lambda = 655; %nm
p.NA = 1.2;
p.pixelsize = 106; %nm
p.sig1=sqrt(2)/2/pi*p.lambda/p.NA/p.pixelsize; %[Zhang 2007]
p.sig2=p.sig1;

pixelizeversion = 0;
vard=zeros(length(l2), length(int1_multi));

for mm=1:length(int1_multi)      
    int_vec=cat(1,int1_multi{mm},int2_multi{mm});   
    % Not integrating out. Cheating...
    % In this case only static situation (int_vec has only one column)
    
    int_vec_static = int_vec/2; % Intensity must be adjusted by 0.5 to reflect the double number of photons compared to the blinking situation. This must be done before adding background. Correcting the final vard by a factor of 2 would not compare the background situation correctly...
    bg_static = p.offset/2; % This is for fair comparison. It is like recording with half integration time. Backgroudn lever is then lower. 
    [vard(:,mm), I3d(:,:,:,mm)]=computeSeparationVariance(x,l1,l2,[p.sig1,p.sig2],int_vec_static, pixelizeversion, bg_static);    
    if length(int1_multi{mm})==1
        % Integrating out 
        [vardintout(:,mm), Iintout(:,:,:,mm)]=computeSeparationVarianceIntOut(x,l1,l2,[p.sig1,p.sig2],int_vec, pixelizeversion, p.offset);
    end
end

% Now with corrected intensity by 0.5 in computaiton of vard we can compare
% vard and vardintout withou any further correction. This takes the effect
% of the background correctly (computation with the same intensity and
% correcting vard by a factor of 2 is not a correct comparison... )