ca
residual = (dpixc-(reshape(res.w*res.h, peval.nx, peval.ny, peval.nt)));
[h, xh] = hist(residual(:)+peval.bg,100);
ixh=[20:65];
[pbest_poiss,perror_poiss,nchi_poiss]=nonlinft('poissoncontinuos' ,xh(ixh),h(ixh)/sum(h(ixh)),sum(h(ixh))./h(ixh),[peval.bg 1.3],[1 1]);
[pbest_gauss,perror_gauss,nchi_gauss]=nonlinft('gaussfunction' ,xh(ixh),h(ixh)/sum(h(ixh)),sum(h(ixh))./h(ixh),[peval.bg 10 1.2],[1 1 1]);
plot(xh,h,'o')
hold on
plot(xh,sum(h)*poissoncontinuos(xh, pbest_poiss),'g')
plot(xh,sum(h)*gaussfunction(xh, pbest_gauss),'r')
grid on
l{1} = 'residuals';
l{2} = 'poisson';
l{3} = 'gauss';
legend(l);
xlabel(residual)
ylabel('count')
xlabel('residual+background')

SaveImageFULL(['residual_hist_nc' num2str(peval.ncomp)])
dipshow(var(residual,[],3))
SaveImageFULL(['residual_var_nc' num2str(peval.ncomp)])
dipshow(mean(residual,3))
SaveImageFULL(['residual_mean_nc' num2str(peval.ncomp)])