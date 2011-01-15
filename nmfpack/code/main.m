function main( what, expnum )
%
% SYNTAX:
% main( what, expnum );
%
% what      - 'learn' or 'show'
% expnum    - experiment number code (specifies data, algorithm, params)
%    
    
global imloadfunc;    
    
% File to store results in
fname = ['../results/experiment-' num2str(expnum) '.mat'];

% Random number generator seed
seed = 0;
randn('seed',seed);
rand('seed',seed);

% Which func loads images? (set to 'imread' to use Matlab's imread.m)
imloadfunc =  'pgma_read'; 

%------------------------------------------------------------------------------
%                            Decode 'expnum'
%------------------------------------------------------------------------------
switch expnum,
  
 % -------------------------- CBCL database --------------------------------   
    
 % Reproduces Lee/Seung Nature-paper (1999) face image results
 case 1,
  dataset = 'cbcl';
  algo = 'nmfdiv';
  rdim = 49;
 
 % Regular NMF but with euclidean objective
 case 2,
  dataset = 'cbcl';
  algo = 'nmfmse';
  rdim = 49;
  
 % Experiments with sparseness constrained NMF
 case 11,
  dataset = 'cbcl';
  algo = 'nmfsc';
  rdim = 49;
  sW = 0.8;
  sH = [];
 case 12,
  dataset = 'cbcl';
  algo = 'nmfsc';
  rdim = 49;
  sW = [];
  sH = 0.8;
 case 13,
  dataset = 'cbcl';
  algo = 'nmfsc';
  rdim = 49;
  sW = 0.2;
  sH = [];
  
  
 % -------------------------- ORL database --------------------------------   
  
 % Regular NMF applied to ORL database => global features
 % (reproduces Li SZ et al results)
 case 101,
  dataset = 'orl';
  algo = 'nmfmse';
  rdim = 25;
  
 % LNMF applied to ORL database => local, almost binary, features
 % (reproduces Li SZ et al results)
 case 111,
  dataset = 'orl';
  algo = 'lnmf';
  rdim = 25;
  alpha = 1; % Note: alpha and beta have no affect on the algorithm,
  beta = 1;  %       they only affect the objective function displayed
  
 % SNMF applied to ORL database, various choices for alpha
 case 121,
  dataset = 'orl';
  algo = 'snmf';
  alpha = 1;
  rdim = 25;
 case 122,
  dataset = 'orl';
  algo = 'snmf';
  alpha = 100;
  rdim = 25;
 case 123,
  dataset = 'orl';
  algo = 'snmf';
  alpha = 10000;
  rdim = 25;
  
 % NMFsc applied to ORL database => local features
 case 131,
  dataset = 'orl';
  algo = 'nmfsc';
  sW = 0.75;
  sH = [];
  rdim = 25;  
 case 133,
  dataset = 'orl';
  algo = 'nmfsc';
  sW = 0.5;
  sH = [];
  rdim = 25;
 case 134,
  dataset = 'orl';
  algo = 'nmfsc';
  sW = 0.6;
  sH = [];
  rdim = 25;
  
  
  
 % ---------------- ON/OFF filtered natural image database ---------------
 
 % Regular NMF, only gives circular 'blobs'
 case 201,
  dataset = 'nat';
  algo = 'nmfmse';
  rdim = 72;
  
 % LNMF, also gives (almost binary) circular 'blobs'
 case 211,
  dataset = 'nat';
  algo = 'lnmf';
  rdim = 72;
  alpha = 1; % Note: alpha and beta have no affect on the algorithm,
  beta = 1;  %       they only affect the objective function displayed
  
 % SNMF, also gives circular 'blobs', with increasing representation
 % error as one increases alpha
 case 221,
  dataset = 'nat';
  algo = 'snmf';
  alpha = 1;
  rdim = 72;
 case 222,
  dataset = 'nat';
  algo = 'snmf';
  alpha = 100;
  rdim = 72;
 case 223,
  dataset = 'nat';
  algo = 'snmf';
  alpha = 10000;
  rdim = 72;
  
 % NMFsc, only one able to produce oriented edge-like features
 case 231,
  dataset = 'nat';
  algo = 'nmfsc';
  sW = [];
  sH = 0.85;
  rdim = 72;
 
    
end



%------------------------------------------------------------------------------
%                             Learn or show
%------------------------------------------------------------------------------

switch what,
  
 % ----------------------- Learn ------------------------
 case 'learn',
  
  % Read in the data and set parameter 
  switch dataset,
   case 'cbcl', V = cbcldata; save('../results/cbcldata.mat','V');
   case 'orl', V = orldata; save('../results/orldata.mat','V');
   case 'nat', V = natdata; save('../results/natdata.mat','V');
  end
  
  % Do not show progress all the time
  showflag = 0;
  
  % Start the algorithm
  switch algo,
   case 'nmfdiv', 
    nmfdiv( V, rdim, fname, showflag );	
   case 'nmfmse', 
    nmfmse( V, rdim, fname, showflag );   
   case 'lnmf',
    lnmf( V, rdim, alpha, beta, fname, showflag );    
   case 'snmf',
    snmf( V, rdim, alpha, fname, showflag );    
   case 'nmfsc',
    nmfsc( V, rdim, sW, sH, fname, showflag ); 
   otherwise,
    error('no such algorithm implemented');	
  end
  
  
 % ----------------------- Show -------------------------
 case 'show',
  
  % Read in the results
  load(fname);

  % Print iteration
  fprintf('Iter: %d\n',iter);
  
  % DISPLAY COLUMNS OF W

  switch dataset,
   case 'cbcl',
    figure(1); 
    normalizecols = 1;
    if normalizecols, 
	Wn = W./(ones(size(W,1),1)*sqrt(sum(W.^2))); 
    else Wn = W;
    end
    visual(Wn,3,7);
   case 'orl',
    figure(1); 
    visual(W,1,ceil(sqrt(size(W,2))),56);
   case 'nat'
    Won = W(1:(end/2),:);
    Woff = W((1:(end/2))+(end/2),:);
    Wdiff = Won-Woff;
    figure(11); 
    visual(Won,2,12);
    figure(12); 
    visual(Woff,2,12);
    figure(13); 
    visual(Wdiff,2,12);
    
  end
  
  % SHOW OBJECTIVE FUNCTION HISTORY
  
  figure(2); 
  if exist('objhistory'),
      if length(objhistory)>1,
	  plot(objhistory(2:end));
      else
	  clf;
      end
  end
  
  % ANALYZE SPARSENESS
  
  vdim = size(W,1);
  samples = size(H,2);
  
  figure(3);
  
  % How much is each unit utilized (in terms of total energy)?
  subplot(3,1,1); bar(sqrt(sum(W.^2)).*sqrt(sum(H'.^2)));
  
  % How sparse are the basis vectors?
  cursW = (sqrt(vdim)-(sum(abs(W))./sqrt(sum(W.^2))))/(sqrt(vdim)-1);
  subplot(3,1,2); bar(cursW);
  
  % How sparse are the coefficients
  cursH = (sqrt(samples)-(sum(abs(H'))./sqrt(sum(H'.^2)))) ...
	  /(sqrt(samples)-1);
  subplot(3,1,3); bar(cursH);
  
  
  % PRINT OUT FINAL APPROXIMATION ERROR

  % Load in the original data
  switch dataset,
   case 'cbcl', load('../results/cbcldata.mat');
   case 'orl', load('../results/orldata.mat');
   case 'nat', load('../results/natdata.mat');
  end

  WH = W*H; 
  Emse = 0.5*sum(sum((V-WH).^2));
  V(find(V<eps))=eps; WH(find(WH<eps))=eps;
  Ediv = sum(sum((V.*log(V./(WH))) - V + WH));
  
  fprintf('---------------------------------\n');
  fprintf('Approximation error:\n');
  fprintf('MSE = %.5f \n',Emse);
  fprintf('Div = %.5f \n',Ediv);
  fprintf('---------------------------------\n');
  
  % SHOW HOW MUCH TIME IT TOOK...
  
  fprintf('Running time (in seconds):\n');
  fprintf('%.1f\n',elapsed);
  fprintf('---------------------------------\n');
  
end


