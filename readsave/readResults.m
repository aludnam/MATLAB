function [res] = readResults(separ, offs, iter, avg, peval)
% [dpixc, dpix_ind, blinkmat, p] = readSimul(separ, offs, iter, avg, p)
% Reads simulations with required 'sep', 'offs', 'iter' and makes
% avgeraging of poisson noise 'avg' times

namedir_data = [peval.prename num2str(100*separ) 'offset_' num2str(offs)];
namefile_data = [namedir_data '-iter_' num2str(iter)];
namefile_res = [namefile_data '_avg' num2str(avg)];
res = load ([peval.path_res namedir_data '/' namefile_res '.mat']);