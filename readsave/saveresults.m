function saveresults(res,peval,p, name_res_file)
% saveresults(res,peval,p, name_res_file)

if ~exist(peval.res_dir,'dir')
    fprintf('Creating directory:\n%s\n', peval.res_dir)
    mkdir (peval.res_dir);
end

if exist('peval.dir_res','var');
    savehere = [peval.dir_res '/' name_res_file];    
else
    savehere = name_res_file;
end

if ~isfield (peval,'fid') peval.fid=1;end
mfprintf(peval.fid, ['Saving data to:\n%s\n'], savehere)
save (savehere,'res', 'peval', 'p');
