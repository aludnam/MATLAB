%                         peval.dir_res_add = ['_iter' num2str(iter)];
% peval.dir_res = [peval.path_res peval.namefile_res];
if ~exist(peval.res_dir,'dir')
    fprintf('Creating directory:\n %s\n', peval.res_dir)
    mkdir (peval.res_dir)
end
if ~exist('p', 'var')
    p=[];
end
saveresults(res, peval, p, [peval.res_path peval.res_dir '/' peval.res_name])
writedata([],[],peval,[peval.res_path peval.res_dir '/' peval.res_name '_param_eval'])
writedata([],[],p,[peval.res_path peval.res_dir '/' peval.res_name '_param_simul'])
if savethis>1 %save intermediate results
%     namedir_tmp = [peval.res_dir '/' 'resetH' num2str(nn) '_iter' num2str(iter)];
    namedir_tmp = [peval.res_dir '/' 'resetH' num2str(nn)];
    mkdir (namedir_tmp)
    cd (namedir_tmp)
    plotreswh2(res, peval, dpixc,p, savethis, 1, 0, 3,1)
    save peval 'peval'
    save res 'res'
    cd (peval.res_path)
end