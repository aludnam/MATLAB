function saveresults(res,peval,p, name_res_file)
% function saveresults(res,peval,p, name_res_file)

% if ~isfield(peval, 'dir_res_add'); peval.dir_res_add=[]; end
% dir_res = [peval.path_res peval.namedir_data peval.dir_res_add];
if ~exist(peval.res_dir,'dir')
    fprintf('Creating directory:\n %s\n', peval.res_dir)
    mkdir (peval.res_dir);
end
% if exist('peval.res_path','var')
%     if ~(strcmp(peval.res_path, peval.data_path)) %not identical
%         fprintf('Creating directory:\n %s\n', peval.dir_res)
%         mkdir (peval.res_dir);
%     end
% end

if exist('peval.dir_res','var');
    savehere = [peval.dir_res '/' name_res_file];    
else
    savehere = name_res_file;
end

if ~isfield (peval,'fid') peval.fid=1;end
mfprintf(peval.fid, ['Saving data to:\n%s\n'], savehere)
save (savehere,'res', 'peval', 'p');
