function [dpixc, dpixc_ind, blinkmat, peval, p] = readSimulSimple(peval)
% [dpixc, dpixc_ind, blinkmat, peval, p] = readSimulSimple(peval)
% Reads simulations with required 'sep', 'offs', 'iter' and makes
% avgeraging of poisson noise 'avg' times
if ~isfield(peval, 'addnamefileres'); peval.addnamefileres=[]; end
if ~isfield(peval, 'catiter_vec'); peval.catiter_vec=[]; end
if ~isfield(peval, 'avg'); peval.avg=1; end

readthis = [peval.data_path peval.data_dir '/' peval.data_name '.mat'];
fprintf ('Reading data: \n')
fprintf ('%s\n', readthis )
load (readthis)

if length(peval.catiter_vec)>1
    fprintf ('Concatenate [')
    fprintf ( '%g ', peval.catiter_vec)
    fprintf ('datasets togoether... \n')
    [dpixc, blinkmat] = catSimul(peval.path_data,  peval.namedir_data, peval.catiter_vec);
    peval.nt = p.Nt*length(peval.catiter_vec);
end

if ~(peval.avg==1) %averaging of Poisson Noise - avg=0 -> noise free
    if peval.avg==0
        fprintf ('Making noise free image \n')
    else
        fprintf ('Making %g averages of each time-slice\n', avg)
    end
    dpixc = makeavg(dpixc_ind, blinkmat, peval.offset, peval.avg);
end