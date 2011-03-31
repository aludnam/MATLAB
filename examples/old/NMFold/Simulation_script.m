peval.path_res = [cd '/'];
peval.path_data = '~/project/data/qdots/S44/';
peval.prename = 'S44_sep_';

peval.sep_vec =[1];
peval.offset_vec = [100];
peval.niter_vec = 1;
peval.avg_num_vec = [1]; %number of averages os samples (reducing Poisson noise)
% avg_num=0 --> noise free image...
peval.ncomp = 2; %number of components without the background component....
peval.catiter_vec=[1]; %cat iterations of simulations - to make data longer

savethis = 1;

for sep_ix=1 : length(peval.sep_vec)                        %separation
    separ = peval.sep_vec(sep_ix);
    for offs_ix = 1: length(peval.offset_vec)               %offset
        offs = peval.offset_vec(offs_ix);
        for iter_ix=1:peval.niter_vec                       %diffrent simulations
            iter = peval.niter_vec(iter_ix);
            for avg_ix=1 : length(peval.avg_num_vec)        %averaging
                avg = peval.avg_num_vec(avg_ix);                
                [dpixc, dpixc_ind, blinkmat, peval, p] = readSimul(separ, offs, iter, avg, peval);
%                 % testing convergence:
%                 dpixc = double(noise(dpixc_ind'*blinkmat(:,372)+p.offset,'poisson'));
                 winit_pix = double(array2im(dpixc_ind));
%                 winit_pix = [];
                hinit = [];
                [res, peval] = separcomp(dpixc, peval, winit_pix, hinit);
                if savethis == 1                    
                    saveresults(res, peval, p, peval.namefile_res)
                    peval.dir_res = [peval.path_res peval.namedir_data];                    
                    writedata([],[],peval,[peval.dir_res '/' peval.namefile_res '_param_eval'])
                    writedata([],[],p,[peval.dir_res '/' peval.namefile_res '_param_simul'])                                        
                end
                
            end %avg
        end %iter
    end %offs
end %separ