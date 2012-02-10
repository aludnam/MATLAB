% clear
savethis = 1;
clear peval

ncomp = [10 20 30]; 


itervec=1; %number of evaluation of each dataset

% set here the path to the data (in this case data must be in mat format):
peval.data_path = '~/project/data/qdots/....';
peval.data_file = 'dpixc';

load([peval.data_path '/' peval.data_file])

% uncomment this for callibration of the data to photon counts:
% bgOffset=700; % this is from the background images
% photonFactor = .2612; %from EM ccd callibration: https://docs.google.com/spreadsheet/ccc?key=0AlBph96P6KPwdEtOcGJWampsaHpGOHBDYUtIZ19LS2c&hl=en_US#gid=0
% dpixc = photonFactor*(dpixc-bgOffset); 
% dpixc=dpixc(:,:,1:300);

% bacground estimation is here:
peval.bg=100; 

[peval.nx, peval.ny, peval.nt]=size(dpixc);
for nc = ncomp;
%     for iterindex=1:niter
    for iterindex=itervec
        if savethis
            [peval.logfile, peval.fid] = initlogfile;
        end
        if sum(dpixc(:)<=0)
            mfprintf(peval.fid,'Clipping negative values in the dpixc!\n')
            dpixc(dpixc<=0)=eps; % to avoid negative values and zeros...
        end

        
        params %reads the parameters        
%         datasource = [peval.home '/' peval.data_path '/' peval.data_dir '/' peval.data_file];
        
%         readdata
        
        % Estimating / subtracting background:
        %         [dpixc, peval]=backgroundestimation(dpixc, peval, p.offset);
        
%         peval.bg=p.offset;
        %         dpixc=bgsubtractbyhand(dpixc,peval);
        peval.ncomp=nc;
        
        dvec=reshape(dpixc,peval.nx*peval.ny,peval.nt);
        
        for indexrestart=0:peval.nrestarts-1
            % Initialization of W:
            %             winit = init_wmap('rand',peval,double(array2im(dpixc_ind)));
            if indexrestart>0
                mfprintf(peval.fid, '\nRestart %g: h restarted and w reused\n',indexrestart);
                winit = init_w_general('res',peval,res.w);
                peval.fix_bg_h=0; mfprintf(peval.fid,'\nBackground component of h will be updated\n')
                
            else            
                winit = init_w_general(peval.init_w_method,peval,[]);
            end
            hinit = init_h_general(peval.init_h_method,peval,[],dpixc);
            
            % Main computation:
            [res, peval]=peval.fun(peval, dvec,winit, hinit);
        end
        
        if savethis
            res.data_file = peval.data_path; % just to keep track where the filel is
            peval = createresname(peval,iterindex);
            if ~exist('p','var')
                p = [];
            end
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
