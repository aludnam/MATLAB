rs = 1;
sizevec = [0 size(dpixc,1), 0 size(dpixc,2)];
Nt = size(dpixc,3);
nx = sizevec(2)-sizevec(1);
ny = sizevec(4)-sizevec(3);
dveccr = squeeze(reshape(dpixc, nx*ny*rs^2, 1, Nt)); % vectors of resized images
numOfIC = 30;
[icasig, A, W] = fastica (dveccr, 'numOfIC', numOfIC, 'g', 'tanh');
sica = size(A,2);
icapix = reshape(A,nx*rs, ny*rs, sica);

fprintf('maxima: \n')
for ii=1:sica
    aim = abs(imresize(icapix(:,:,ii),1/rs));
    m(ii) = max(max(aim));
    [ym(ii), xm(ii)] = find(m(ii)==aim);
    [ycog(ii), xcog(ii)] = cog(imresize(icapix(:,:,ii),1/rs)); %center of gravity
    [xfit(ii), yfit(ii), sig(ii)] = fitgauss2d(imresize(icapix(:,:,ii),1/rs));
    fprintf('%g . maximum: %g %g \tcog: %g %g \tfit: %g %g \tsig: %g\n',ii , [xm(ii), ym(ii)], [xcog(ii) ycog(ii)], [xfit(ii), yfit(ii), sig(ii)] );
end

% goodones = find(and(sig<s+1, sig>s-1));
% resfit = pixelize(round([xfit; yfit]'), ones(length(yfit),1), [0, nx, 0, ny],nx, ny,[], 0);
% resfit_good = pixelize(round([xfit(goodones); yfit(goodones)]'), ones(length(goodones),1), [0, nx, 0, ny],nx, ny,[], 0);
% % rescog = pixelize(round([xcog; ycog]'), ones(length(xcog),1), [0, nx, 0, ny],nx, ny,[], 0);
% 
% % ims((res(12:32, 11:31)))
% res = resfit_good;
% 
% ims(res);
ims(sum(dpixc,3));

imstiled(imresize(icapix,1),[],0)
figure;
ims(sum(dpixc,3),'gray');
hold on
plotData([xfit; yfit]'*rs,[0 nx 0 ny]*rs,'xb',20);
hold off
