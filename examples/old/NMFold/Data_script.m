peval.path_res = [cd '/'];
peval.path_data = '~/project/data/qdots/S44/';

peval.ncomp = 2; %number of components without the background component....

savethis = 1;
winit_pix = [];
hinit = [];
[res, peval] = separcomp(dpixc, peval, winit_pix, hinit);
if savethis == 1
    saveresults(res, peval, p, peval.namefile_res)
    peval.dir_res = [peval.path_res peval.namedir_data];
    writedata([],[],peval,[peval.dir_res '/' peval.namefile_res '_param_eval'])
    writedata([],[],p,[peval.dir_res '/' peval.namefile_res '_param_simul'])
end
