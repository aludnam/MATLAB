function saveresults2(res,peval)
fprintf('saving data \n');
dir_res = [peval.path_res peval.name_res_dir];
if ~(strcmp(peval.path_res, peval.path_data)) %not identical
    mkdir (dir_res);
end
save ([dir_res '/' peval.namefile_res],'res', 'peval')