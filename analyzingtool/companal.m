function [xm, ym, xfit,yfit] = companal(icapix, rs)

sica = size(icapix,3);
fprintf('maxima: \n')
for ii=1:sica
    aim = abs(imresize(icapix(:,:,ii),1/rs));
    m(ii) = max(max(aim));
    [ym(ii), xm(ii)] = find(m(ii)==aim);
    [ycog(ii), xcog(ii)] = cog(imresize(icapix(:,:,ii),1/rs)); %center of gravity
    [xfit(ii), yfit(ii), sig(ii)] = fitgauss2d(imresize(icapix(:,:,ii),1/rs));
    fprintf('%g . maximum: %g %g \tcog: %g %g \tfit: %g %g \tsig: %g\n',ii , [xm(ii), ym(ii)], [xcog(ii) ycog(ii)], [xfit(ii), yfit(ii), sig(ii)] );
end
