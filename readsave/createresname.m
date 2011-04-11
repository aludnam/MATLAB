function peval = createresname(peval, app)
% peval = createresname(peval, app)
    
    
res_nameappendix = ['_nc' num2str(peval.ncomp)];
if isfield(peval, 'data_dir')
    peval.res_dir = [peval.data_dir res_nameappendix];
else
    if isfield(peval, 'fun')
        res_nameappendix = ['_' func2str(peval.fun) res_nameappendix];
    end
    peval.res_dir = ['results' res_nameappendix];
end

peval.res_name = ['results_iter' num2str(app)];