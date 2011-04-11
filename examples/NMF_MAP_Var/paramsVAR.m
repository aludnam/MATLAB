% params
mfprintf(peval.fid, 'Parameters read from:\n%s\n',mfilename('fullpath'));

peval.fun = @updates_variational;
% peval.fun = @variationalupdates;
peval.home = getenv('HOME');

% Reading data
peval.data_path = 'project/data/qdots/D3ptBMrand';
peval.data_dir = data_core;
peval.data_file = [data_core '_iter1.mat'];

peval.res_path = [cd '/'];
peval.maxiter = 500;        %maximum number of iteration
peval.init_w_method = 'rand'; % initilaization method for W
peval.init_h_method = 'rand'; % initilaization method for H

peval.bgcomp = 1;   % last component as a background
peval.fix_bg_w = 1; % fixes last component w in the evaluation
peval.fix_bg_a = 1; % fixes last component a in the evaluation

% Prior parameters of the Gamma distsribution
peval.alpha = 1; 
peval.beta = 1/1000;

peval.showprogress = 1; %plots the progress of the evaluation

