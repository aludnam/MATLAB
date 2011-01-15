function saveresults(res,peval,p, name_res_file)
% function saveresults(res,peval,p, name_res_file)

% if ~isfield(peval, 'dir_res_add'); peval.dir_res_add=[]; end
% dir_res = [peval.path_res peval.namedir_data peval.dir_res_add];
if exist('peval.path_res','var')
    if ~(strcmp(peval.path_res, peval.path_data)) %not identical
        fprintf('Creating directory:\n %s\n', peval.dir_res)
        mkdir (peval.dir_res);
    end
end

if exist('peval.dir_res','var');
    savehere = [peval.dir_res '/' name_res_file];    
else
    savehere = name_res_file;
end
fprintf('Saving data to: \n %s \n', savehere);
save (savehere,'res', 'peval', 'p');
