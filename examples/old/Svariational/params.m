% params
mfprintf(peval.fid, 'Parameters read from:\n%s\n',[cd '/params.m'])

peval.home = getenv('HOME');

% Reading data
data_N = 2;
data_offset = 100;
data_iter = 1;

peval.data_path = 'project/data/qdots/D2ptBMrand';
peval.data_dir = data_core;
peval.data_file = [data_core '.mat'];

peval.res_path = [cd '/'];
peval.maxiter = 250;        %maximum number of iteration
peval.initw1_method = 'rand'; % initilaization method for W

% Prior parameters of the Gamma distsribution
peval.alpha = 1; 
peval.beta = 1/1000;

peval.bg_clip = 'yes'; % clips negative and zero values when subtracting background
peval.addbgcomp = 1; % adds a background component
peval.showprogress = 1; %plots the progress of the evaluation
