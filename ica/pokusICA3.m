% tstart=tic;
cd ~/project/MATLAB/PatternAnalysis/
%N = 2;
Nt = 10000;
rs = 0.5; %resizinf faction
sizevec = [0 128, 0 128];
nx = sizevec(2)-sizevec(1);
ny = sizevec(4)-sizevec(3);
% N = 20; d = generateNoise(N, sizevec);
%d = [2 20; 2 2; 2 5; 20 2; 20 5; 20 8];
% d = [10 8; 11, 8; 12, 8; 13, 8; 14, 8; 10 10; 11, 10; 12, 10; 13, 10; 14, 10];

% d(1:15,1)=51:65;
% d(1:15,2)=50;
% d(16:30,1)=51:65;
% d(16:30,2)=53;
% 
% 
% d(31:61,1) = 3*sin(0.5*[1:31]) + 50;
% d(31:61,2) = 61:91;
% 
% d(62:92,1) = 3*sin(0.5*[1:31]) + 53;
% d(62:92,2) = 61:91;
% 
% d(62:92,1) = 3*sin(0.5*[1:31]) + 56;
% d(62:92,2) = 61:91;
% 
% d(93:123,1) = 5*sin(0.5*[1:31]) +80;
% d(93:123,2) = 5*cos(0.5*[1:31]) +60;
% 
% d(124:144,1) = 2*sin(0.5*[1:21]) +80;
% d(124:144,2) = 2*cos(0.5*[1:21]) +60;

d=[64 64; 63 64];

N = size(d,1);
% blinkmat = rand(N, Nt)>0.5;
blinkmat = rand(N, Nt);
dpix_orig = pixelize(d, ones(size(d,1),1), sizevec, nx, ny ,[], 0);

lambda = 250; %nm
pixelsize = 106; %nm
lambdapix = lambda/pixelsize; %pixels
NA = 1.2;
s = 8.8*NA/lambdapix; %gaussian approx of airy
% s = 1/rs*s; %for before resizing...

[X, Y] = meshgrid(-20:20);
f = exp( -(X.^2+Y.^2)/(2*s^2) );
f = f / sum(f(:));
sf = size(f)-1;

for ii=1 : N
    dpix_ind(:,:,ii) = pixelize(d(ii,:), 1, sizevec, nx, ny, [],0);
    dpixc_ind(:,:,ii) = conv2(dpix_ind(:,:,ii),f,'same'); %convolution of individual points
    
end

dvecc_ind = squeeze(reshape(dpixc_ind, nx*ny, 1, N)); % each indiv image to vector

dvecc = dvecc_ind*blinkmat; % set of Nt vectors
dpixc = imresize(reshape(dvecc, nx, ny, Nt), rs); % resized images from vectors dvecc


% toc(tstart)
cd ~/project/MATLAB/FastICA_25/
numOfIC = N;
[icasig, A, W] = FASTICA (dvecc', 'numOfIC', numOfIC, 'g', 'tanh');
% [icasig] = FASTICA (dvec', 'numOfIC', numOfIC , 'g', 'pow3');

sica = size(icasig,1);
% icapix = shiftdim(reshape(icasig,sica, nx+sf(1), ny+sf(2)),1);
icapix = shiftdim(reshape(icasig,sica, nx*rs, ny*rs),1);

for ii=1:sica
    aim = abs(imresize(icapix(:,:,ii),1/rs));
    m(ii) = max(max(aim));
    [ym(ii), xm(ii)] = find(m(ii)==aim);
end
ym = ny - ym; %just corect

cd ~/project/MATLAB/PatternAnalysis/
% res = pixelize([xm;ym]', ones(length(xm)), [0, nx+sf(1), 0, ny+sf(2)],nx+sf(1) , ny+sf(2) ,[], 0);
res = pixelize([xm;ym]', ones(length(xm),1), [0, nx, 0, ny],nx, ny,[], 0);
cd ~/project/MATLAB/


% ims((res(12:32, 11:31)))
ims(res)
ims(sum(dpix,3))
ims(sum(dpixc,3))