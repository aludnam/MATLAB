function bicplot(ncomp_vec,bicp,loglik,penalty)
% function bicplot(ncomp_vec,bicp,loglik,penalty)
figure
[AX,H1,H2] = plotyy(repmat(ncomp_vec',1,2),[bicp',loglik'], ncomp_vec, penalty);
set(get(AX(1),'Ylabel'),'String','BIC / log(likelihood)')
set(get(AX(2),'Ylabel'),'String','penalty')
grid on
xlabel('number of components')
l={'BIC','log(likel)','penalty'};
legend(l)