function saveresults_data(res,peval, name_res_file)
fprintf('saving data \n');
dir_res = [peval.path_res];
if ~(strcmp(peval.path_res, peval.path_data)) %not identical
    mkdir (dir_res);
end
save ([dir_res '/' name_res_file],'res', 'peval', 'p')