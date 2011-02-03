
%                         peval.dir_res_add = ['_iter' num2str(iter)];
peval.dir_res = [peval.path_res peval.namefile_res];
if ~exist(peval.dir_res,'dir')
    fprintf('Creating directory:\n %s\n', peval.dir_res)
    mkdir (peval.dir_res)
end
saveresults(res, peval, p, [peval.dir_res '/' peval.namefile_res])
writedata([],[],peval,[peval.dir_res '/' peval.namefile_res '_param_eval'])
writedata([],[],p,[peval.dir_res '/' peval.namefile_res '_param_simul'])
if savethis>1 %save intermediate results
    namedir_tmp = [peval.dir_res '/' 'resetH' num2str(nn) '_iter' num2str(iter)];
    mkdir (namedir_tmp)
    cd (namedir_tmp)
    plotreswh2(res, peval, dpixc,p, savethis, 1, 0, 3,1)
    save peval 'peval'
    cd (peval.path_res)
end