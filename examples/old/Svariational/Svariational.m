clear
savethis = 1;
% delta = [0.1:0.1:1,1.2];
delta = [0.2 .6];
niter=2; %number of evaluation of each dataset
for deltaindex=1:length(delta)
    for iterindex=1:niter
        ca;
        data_core = ['D_N2_delta' num2str(delta(deltaindex)*100)];        
        if savethis
            peval.logfile = 'evaluation.log';
            fidlog = logfileinit(peval.logfile);
            peval.fid = [1 fidlog];
        end
        
        %reads the parameters
        paramsvariational
        datasource = [peval.home '/' peval.data_path '/' peval.data_dir '/' peval.data_file];
        % reads the data
        readdata
                        
        % Estimating / subtracting background:
%         [dpixc, peval]=backgroundestimation(dpixc, peval, p.offset);        
        peval.bg=p.offset;                
        
%         peval.ncomp=p.N;
        peval.ncomp=3;
        
        % Initialization of W:
        winit = init_w(peval.initw1_method,peval,double(array2im(dpixc_ind)));
        
        % need fit to blinkmat (1xpeval.ncomp)
        alph=repmat(peval.alpha,1,peval.ncomp);
        beta=repmat(peval.beta,1,peval.ncomp);
        
        % initialization of a_k, b_k
        [ainit, binit] = init_ab(dpixc, alph, beta, peval);
       
        [res, peval]=variationalupdates(peval, dvec,winit, alph,beta, ainit, binit);        
        
        if savethis
            peval = createresname(peval,iterindex);            
            saveresults(res, peval, p, [peval.res_path peval.res_dir '/' peval.res_name])            
            fclose(peval.fid(2));
        end        
    end
    if savethis        
        % Saving parametes and logfile into the directory with resuls
        saveparameters(peval,p)
    end
end

% plottingmulti
% plotscattehvsblinkmat