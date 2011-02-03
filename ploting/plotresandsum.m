function plotresandsum(dpixc,w,peval,savethis,p)
% plotresandsum(dpixc,w,peval,savethis,p)

dipshow(mean(dpixc,3))
hold on
[x_mu, y_mu sig] = localizew(w,peval);
scatter(x_mu-1, y_mu-1,'r')
if exist('p','var')
    scatter(p.x_vec, p.y_vec,'xb')
end
drawpixels([peval.nx, peval.ny],':r',-0.5)

if savethis
    SaveImageFULL([peval.res_dir '/localizedw'])
end