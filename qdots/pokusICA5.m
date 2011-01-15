clear
tstart=tic;
addpath '~/project/MATLAB/PatternAnalysis/'
addpath '~/project/MATLAB/qdots/'
addpath '~/project/MATLAB/qdots/FastICA_25/'
%N = 2;
Nt = 500;
rs = 0.25; %resizinf faction
sizevec = [0 32, 0 32];
nx = sizevec(2)-sizevec(1);
ny = sizevec(4)-sizevec(3);
% N = 30; d = generateNoise(N, [sizevec(1)+10, sizevec(2)-10, sizevec(3)+10, sizevec(4)-10]);
%d = [2 20; 2 2; 2 5; 20 2; 20 5; 20 8];
% d = [10 8; 11, 8; 12, 8; 13, 8; 14, 8; 10 10; 11, 10; 12, 10; 13, 10; 14, 10];
% d(:,1)=d(:,1)+20; 
% d(:,2)=d(:,2)+20;
% d(1:15,1)=51:65;
% d(1:15,2)=50;
% d(16:30,1)=51:65;
% d(16:30,2)=53;


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
%

% % x = [1:20];
% % sx = size(x,2);
% % r1 = 10;
% % c = [80, 60]; 
% % d(1:sx,1) = r1 * sin((2*pi/sx) * x) + c(1);
% % d(1:sx,2) = r1 * cos((2*pi/sx) * x) + c(2);
% % 
% % x2 = [1:10];
% % sx2 = size(x2,2);
% % r2 = 5;
% % d(sx+1: sx+sx2,1) = r2 * sin((2*pi/sx2) * x2) + c(1);
% % d(sx+1: sx+sx2,2) = r2 * cos((2*pi/sx2) * x2) + c(2);
 
% d = [80, 60];
% d=[64 64; 63 64; 62 64; 61 64];
%  d=[64 64; 65 64; 64 65; 15 20];
center = ceil([nx, ny]/2);
d = [center; center + [1 0]; center + [2 0]];
% d=[64 64; 65 64; 66 64];

N = size(d,1);

% blinkmat = rand(N, Nt)>0.5;
blinkmat = rand(N, Nt);
dpix_orig = pixelize(d, ones(size(d,1),1), sizevec, nx, ny ,[], 0);

lambda = 400; %nm
pixelsize = 106*rs; %in nm after resizing... 
lambdapix = lambda/pixelsize; %pixels
NA = 1.0;
s = 1.4/(2*pi)*(lambdapix/NA); %gaussian approx of airy
% s = 1/rs*s; %for before resizing...

[X, Y] = meshgrid(-20:20);
psf = exp( -(X.^2+Y.^2)/(2*s^2) );
psf = psf / sum(psf(:));
spsf = size(psf)-1;

maxphot = 100; % maximal expexted number of photons in one pixel 
offset = 0.1; % general offset as a fraction of maximum

% pixelization of each point
for ii=1 : N
    dpix_ind(:,:,ii) = pixelize(d(ii,:), 1, sizevec, nx, ny, [],0);
    dpixc_ind(:,:,ii) = conv2(dpix_ind(:,:,ii),psf,'same'); %convolution of individual points    
end

% pmid = zeros(nx, ny);
% pmid (floor(nx/2), floor (ny/2)) = 1;
% pmid = imresize(conv2(pmid,f,'same'),rs);
% pmidvec = reshape(pmid, nx*ny*rs^2,1);
% guess = repmat (pmidvec,1,N);

dvecc_ind = squeeze(reshape(dpixc_ind, nx*ny, 1, N)); % each indiv image to vector

dpix_ind_small = imresize(dpixc_ind,rs);
dvecc_ind_small = squeeze(reshape(dpix_ind_small, nx*ny*rs^2, 1, N)); % each indiv image to vector

dvecc_nonoise = (dvecc_ind*blinkmat); % set of Nt vectors
dpixc_nonoise = imresize(reshape(dvecc_nonoise, nx, ny, Nt), rs);
dpixc_nonoise_dip = dip_image(dpixc_nonoise);

maxdvec = max(dpixc_nonoise(:));
offset_abs = offset/(1-offset)*maxdvec;
dpixc_dip = noise((dpixc_nonoise_dip+offset_abs)/(maxdvec+offset_abs)*maxphot,'poisson');
dpixc = double(dpixc_dip);
% dpixc = imnoise(uint16((dpixc_nonoise + offset_abs)/(maxdvec + offset_abs)*maxphot),'poisson');
% dpixc = imnoise((dpixc_nonoise + offset_abs)/(maxdvec + offset_abs)*maxphot,'poisson'); 
% dpixc = imresize(reshape(dvecc, nx, ny, Nt), rs); % resized images from vectors dvecc
dveccr = double(squeeze(reshape(dpixc, nx*ny*rs^2, 1, Nt))); % vectors of resized images

toc(tstart)

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

fprintf('maxima: \n')
for ii=1:sica
    aim = abs(imresize(icapix(:,:,ii),1/rs));
    m(ii) = max(max(aim));
    [ym(ii), xm(ii)] = find(m(ii)==aim);
    [ycog(ii), xcog(ii)] = cog(imresize(icapix(:,:,ii),1/rs)); %center of gravity
    [xfit(ii), yfit(ii), sig(ii)] = fitgauss2d(imresize(icapix(:,:,ii),1/rs));
    fprintf('%g . maximum: %g %g \tcog: %g %g \tfit: %g %g \tsig: %g\n',ii , [xm(ii), ym(ii)], [xcog(ii) ycog(ii)], [xfit(ii), yfit(ii), sig(ii)] );
end

goodones = find(and(sig<s+1, sig>s-1));
resfit = pixelize(round([xfit; yfit]'), ones(length(yfit),1), [0, nx, 0, ny],nx, ny,[], 0);
resfit_good = pixelize(round([xfit(goodones); yfit(goodones)]'), ones(length(goodones),1), [0, nx, 0, ny],nx, ny,[], 0);
% rescog = pixelize(round([xcog; ycog]'), ones(length(xcog),1), [0, nx, 0, ny],nx, ny,[], 0);

% ims((res(12:32, 11:31)))
res = resfit_good;

% ims(res);
% ims(sum(dpix_ind,3));
ims(sum(dpixc,3));
ims(var(dpixc,0,3));

% ims(imresize(sum(dpix_ind,3),rs));
% ims(imresize(res,rs));
ims(sum(dpix_ind,3)-res);
imstiled(imresize(icapix,1),[],0)
imstiled(imresize(icapix,1/rs),[],0)
figure;
plotData([xfit; yfit]',[0 nx 0 ny],'xb',20)
hold on; 
plotData(d,[0 nx 0 ny],'r',20)
hold off;
ims(sum(dpixc,3),'gray');
hold on
plotData([xfit; yfit]'*rs,[0 nx 0 ny]*rs,'xb',20);
hold off
toc(tstart)
% ims(imresize(sum(dpix_ind,3),rs)-imresize(res,rs))

% for ii=1 :sica
%     rf(:,:,ii) = pixelize(round([xfit(ii); yfit(ii)]'), 1, [0, nx, 0, ny],nx, ny,[], 0);
%     rf(:,:,ii) = conv2(rf(:,:,ii),f,'same');
% end
% rfs = squeeze(reshape(imresize(rf,rs),nx*ny*rs^2,1,sica));
% [icasig, A, W] = fastica (dveccr, 'numOfIC', numOfIC, 'g', 'tanh','initGuess', rfs);