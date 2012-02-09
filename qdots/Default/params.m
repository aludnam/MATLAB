 
% params
mfprintf(peval.fid, 'Parameters read from:\n%s\n',mfilename('fullpath'));

peval.fun = @updates_nmfclassic;
% peval.fun = @variationalupdates;
peval.home = getenv('HOME');

% Reading data
% peval.data_path = 'project/data/qdots/D3ptBMrand';
% peval.data_dir = data_core;
% peval.data_file = [data_core '_iter1.mat'];

peval.res_path = [cd '/'];
peval.maxiter = 1000;        %maximum number of iteration
% peval.init_w_method = 'image_repmat'; % initilaization method for W
% peval.init_w_method = 'res';
peval.init_w_method = 'rand';
peval.init_h_method = 'rand'; % initilaization method for H

peval.nrestarts = 3; % number of partial restarts (restarting H but keeping W) 

peval.bgcomp = 1;   % last component as a background
peval.fix_bg_w = 1; % fixes last component w in the evaluation
peval.fix_bg_h = 0; % fixes last component a in the evaluation

peval.showprogress = 1; %plots the progress of the evaluation

