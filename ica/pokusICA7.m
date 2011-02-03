clear
tstart=tic;
addpath '~/project/MATLAB/PatternAnalysis/'
addpath '~/project/MATLAB/qdots/'
addpath '~/project/MATLAB/qdots/FastICA_25/'
Nt = 500;
rs = 1; %resizing fraction
sizevec = [0 32, 0 32];
nx = sizevec(2)-sizevec(1);
ny = sizevec(4)-sizevec(3);

lambda = 400; %nm
pixelsize = 106*rs; %in nm after resizing...
lambdapix = lambda/pixelsize; %pixels
NA = 1.6;
cheatfactor = 3; %NA is not really that good
s = cheatfactor*1.4/(2*pi)*(lambdapix/NA); %gaussian approx of airy
% s = 1/rs*s; %for before resizing...

[X, Y] = meshgrid(-20:20);
psf = exp( -(X.^2+Y.^2)/(2*s^2) );
psf = psf / sum(psf(:));

maxphot = 100; % maximal expexted number of photons in one pixel 
offset = 0.1; % general offset as a fraction of maximum

[dpixc, dveccr, N] = generatedata(d, sizevec, psf, maxphot, offset, Nt, rs);

numOfIC = N;

% [p, Ap, Wp] = fastica (dveccr, 'numOfIC', numOfIC, 'g', 'tanh', 'only', 'pca');
% % [icasig, A, W] = fastica (dveccr, 'numOfIC', numOfIC, 'g', 'tanh');
% [icasig, A, W] = fastica (dveccr, 'numOfIC', numOfIC, 'g', 'tanh','approach', 'symm')
% [icasig] = FASTICA (dvec', 'numOfIC', numOfIC , 'g', 'pow3');
% guess = repmat(11);
%  [icasig, A, W] = fastica (dveccr, 'numOfIC', numOfIC, 'g', 'tanh','initGuess', guess);
[icasig, A, W] = fastica (dveccr, 'numOfIC', numOfIC, 'g', 'tanh');
sica = size(A,2);
% icapix = shiftdim(reshape(icasig,sica, nx+sf(1), ny+sf(2)),1);
icapix = reshape(A,nx*rs, ny*rs, sica);

%%% nmf
[w,h]=nmf(dveccr',numOfIC,1);
icapix=shiftdim(reshape(h,numOfIC,nx*rs,ny*rs), 1);

[xm, ym, xfit,yfit] = companal(icapix, rs);

ims(sum(dpixc,3));
ims(var(dpixc,0,3));

imstiled(imresize(icapix,1),[],0)
imstiled(imresize(icapix,1/rs),[],0)
sizedata = 40;
figure;
plotData([xfit; yfit]',[0 nx 0 ny],'xb',sizedata)
hold on; 
plotData(d,[0 nx 0 ny],'r',sizedata)
hold off;
ims(sum(dpixc,3),'gray');
hold on
plotData([xfit; yfit]'*rs,[0 nx 0 ny]*rs,'xb',sizedata);
hold off
toc(tstart)
% ims(imresize(sum(dpix_ind,3),rs)-imresize(res,rs))

% for ii=1 :sica
%     rf(:,:,ii) = pixelize(round([xfit(ii); yfit(ii)]'), 1, [0, nx, 0, ny],nx, ny,[], 0);
%     rf(:,:,ii) = conv2(rf(:,:,ii),f,'same');
% end
% rfs = squeeze(reshape(imresize(rf,rs),nx*ny*rs^2,1,sica));
% [icasig, A, W] = fastica (dveccr, 'numOfIC', numOfIC, 'g', 'tanh','initGuess', rfs);