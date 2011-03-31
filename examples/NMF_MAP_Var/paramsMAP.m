% params
mfprintf(peval.fid, 'Parameters read from:\n%s\n',[cd '/params.m'])

peval.fun = @updates_map;
% peval.fun = @variationalupdates;
peval.home = getenv('HOME');

% Reading data
peval.data_path = 'project/data/qdots/D2ptBMrand';
peval.data_dir = data_core;
% peval.data_file = [data_core '_iter1.mat'];
peval.data_file = [data_core];

peval.res_path = [cd '/'];
peval.maxiter = 500;        %maximum number of iteration
peval.initw1_method = 'rand'; % initilaization method for W

peval.bg_clip = 'yes'; % clips negative and zero values when subtracting background

peval.showprogress = 1; %plots the progress of the evaluation

peval.ddterm = eps;
peval.maxiter = 500;
peval.maxh = 10;
peval.maxw = 10;

peval.init_w_method = 'res';
peval.init_h_method = 'rand';
peval.bgcomp=0;

% % % MAP fitting part parameters
peval.Witer=3;          %number of updates of W
peval.Hiter=20;         %number of updates of H
peval.nIterAlter = 20;  %number of cycles of updates W and H

peval.Wupdate='loglikGaPExpHfix';
peval.Hupdate='loglikGaPExpWfix';
% Prior parameters of the Gamma distsribution
% beta parameter is reciprocal to the one in the variational updates...
peval.alpha = 1;        %exponentional prior on H
peval.beta=1000;

peval.sigmaPSF=1;



