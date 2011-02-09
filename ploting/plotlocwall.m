function plotlocwall(x_mu,y_mu,peval,datalocation,savethis)
% plotlocwall(x_mu,y_mu,peval,datalocation,savethis)
load (datalocation)
% load (['~/' peval.data_path '/' peval.data_dir '/' peval.data_file])
dipshow(mean(dpixc,3))
hold on
scatter(x_mu(:,1)-1, y_mu(:,1)-1, 'dg')
scatter(x_mu(:,2)-1, y_mu(:,2)-1, 'dg')
scatter(p.x_vec,p.y_vec,'xb')
hold off
if savethis
    SaveImageFULL([peval.res_dir '/localizedwall'])
end