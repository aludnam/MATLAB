cd ~/project/MATLAB/PatternAnalysis/
%N = 2;
Nt = 5;
sizevec = [0 128, 0 128];
nx = sizevec(2)-sizevec(1);
ny = sizevec(4)-sizevec(3);
% N = 20; d = generateNoise(N, sizevec);
%d = [2 20; 2 2; 2 5; 20 2; 20 5; 20 8];
% d = [10 8; 11, 8; 12, 8; 13, 8; 14, 8; 10 10; 11, 10; 12, 10; 13, 10; 14, 10];

% % % d(1:15,1)=51:65;
% % % d(1:15,2)=50;
% % % d(16:30,1)=51:65;
% % % d(16:30,2)=53;
% % % 
% % % 
% % % d(31:61,1) = 3*sin(0.5*[1:31]) + 50;
% % % d(31:61,2) = 61:91;
% % % 
% % % d(62:92,1) = 3*sin(0.5*[1:31]) + 53;
% % % d(62:92,2) = 61:91;
% % % 
% % % d(62:92,1) = 3*sin(0.5*[1:31]) + 56;
% % % d(62:92,2) = 61:91;
% % % 
% % % d(93:123,1) = 5*sin(0.5*[1:31]) +80;
% % % d(93:123,2) = 5*cos(0.5*[1:31]) +60;
% % % 
% % % d(124:144,1) = 2*sin(0.5*[1:21]) +80;
% % % d(124:144,2) = 2*cos(0.5*[1:21]) +60;

d=[64 64; 62 62];

N = size(d,1);
blinkmat = rand(N, Nt)>0.5;
% blinkmat = rand(N, Nt);
dpix_orig = pixelize(d, ones(size(d,1),1), sizevec, nx, ny ,[], 0);

s = 3;
[X, Y] = meshgrid([-10:10]);
f = exp( -(X.^2+Y.^2)/(2*s^2) );
f = f / sum(f(:));
sf = size(f)-1;

for ii=1:Nt
    dblink{ii} = d(blinkmat(:,ii),:);
    if isempty (dblink{ii})
        dpix(:,:,ii) = zeros(nx, ny);
    else 
        intens = ones(size(dblink{ii}, 1), 1);
        dpix(:,:,ii) = pixelize(dblink{ii}, intens, sizevec, nx, ny ,[], 0);
    end
    dpixc(:,:,ii) = conv2(dpix(:,:,ii),f,'same');
end

% dvec = squeeze(reshape(dpixc, (nx+sf(1))*(ny+sf(2)), 1, Nt));
dvec = squeeze(reshape(dpixc, (nx)*(ny), 1, Nt));

cd ~/project/MATLAB/FastICA_25/
numOfIC = N;
[icasig] = FASTICA (dvec', 'numOfIC', numOfIC, 'g', 'tanh');
% [icasig] = FASTICA (dvec', 'numOfIC', numOfIC , 'g', 'pow3');

sica = size(icasig,1);
% icapix = shiftdim(reshape(icasig,sica, nx+sf(1), ny+sf(2)),1);
icapix = shiftdim(reshape(icasig,sica, nx, ny),1);

for ii=1:sica
    aim = abs(icapix(:,:,ii));
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