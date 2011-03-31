function saveparameters(peval,p)
% saveparameters(peval)
% saves parameters and logfile to the directory with results

copyfile (peval.logfile, peval.res_dir);
if isfield(peval,'fun') && isa(peval.fun,'function_handle')
        peval.fun=func2str(peval.fun); %fprintf cant write handles....
end
writedata([],[],peval,[peval.res_path peval.res_dir '/param_eval'])
writedata([],[],p,[peval.res_path peval.res_dir '/param_simul'])