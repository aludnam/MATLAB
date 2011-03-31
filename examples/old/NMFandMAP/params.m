% params
fprintf('Parameters read from:\n%s\n',[cd '/params.m'])
peval.res_path = [cd '/'];
% peval.sep_vec =[1];
peval.N = 2;
peval.offset = 100;
% (reducing Poisson noise)
% avg_num=0 --> noise free image...
peval.ddterm = eps;
peval.maxiter = 500;
peval.maxh = 10;
peval.maxw = 10;
peval.method = 'nmf_classic';
peval.nrestartH = 2; %number of restarts of H

peval.initw1_method = 'image_repmat'; %first intialize with mean(dpixc,3)
peval.initw_method = 'res';  % then initialize with (res.w) from before

peval.w_fixvec = []; %figing only added backround component
peval.h_fixvec = [];

% % % MAP fitting part parameters
peval.Witer=3;          %number of updates of W
peval.Hiter=20;         %number of updates of H
peval.nIterAlter = 5;  %number of cycles of updates W and H

peval.Wupdate='loglikGaPExpHfix';
peval.Hupdate='loglikGaPExpWfix';
peval.sigmaPSF = p.sigmapix; %sigma of the PSF gauss approx (in pixels)
% peval.alpha = 1;        %exponentional prior on H
% peval.beta=1000;

%optimization parameters
optionsW = zeros(1,18);
optionsW(1)=1;                  %to display error values
optionsW(7)=1;
optionsW(9)=0;                   %to check gradient
optionsW(14)=peval.Witer;         %maximum number of iterations

optionsH=optionsW;
optionsH(14)=peval.Hiter;
