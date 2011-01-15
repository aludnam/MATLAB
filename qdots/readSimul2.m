function [dpixc, dpixc_ind, blinkmat, peval, p] = readSimul(separ, offs, iter, avg, peval)
% [dpixc, dpix_ind, blinkmat, p] = readSimul(separ, offs, iter, avg, p)
% Reads simulations with required 'sep', 'offs', 'iter' and makes
% avgeraging of poisson noise 'avg' times

peval.namedir_data = [peval.prename num2str(100*separ) 'offset_' num2str(offs)];
peval.namefile_data = [peval.namedir_data '-iter_' num2str(iter)];
peval.namefile_res = [peval.namefile_data '_avg' num2str(avg)];
load ([peval.path_data peval.namedir_data '/' peval.namefile_data '.mat'])

if length(peval.catiter_vec)>1
    fprintf ('Concatenate [')
    fprintf ( '%g ', peval.catiter_vec)
    fprintf ('datasets togoether... \n')
    [dpixc, blinkmat] = catSimul(peval.path_data,  peval.namedir_data, peval.catiter_vec);
    p.Nt = p.Nt*length(peval.catiter_vec);
end

if ~(avg==1) %averaging of Poisson Noise - avg=0 -> noise free
    if avg==0
        fprintf ('Making noise free image \n')
    else
        fprintf ('Making %g averages of each time-slice\n', avg)
    end
    dpixc = makeavg(dpixc_ind, blinkmat, offs, avg);
end