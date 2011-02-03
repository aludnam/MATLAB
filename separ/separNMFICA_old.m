% evaluate simulated data - ICA and NMF
function separNMFICA(sep0, offset0, path_data, path_res, prename, niter, savethis, initval, sep_how)

if ~exist('initval', 'var')
    initval = 0;
end

if ~exist('sep_how', 'var')
    sep_how = 'in';
end
ncomp = 3; %number of components to be separated + background

for rr = 1: length(offset0)
    for ll=1 : length(sep0)
        namedir = [prename num2str(100*sep0(ll)) 'offset_' num2str(offset0(rr))];
        cd ([path_data namedir])
        for mm=1:niter
            
            namefile = [namedir '-iter_' num2str(mm)];
            load ([namefile '.mat'])
            
            if sum(sep_how == 'i')>0 %ICA
                [icasig{mm}, A{mm}, W{mm}] = fastica (dveccr, 'numOfIC', ncomp, 'g', 'tanh');
                icapixICA{mm} = reshape(A{mm},p.nx, p.ny, ncomp);
            end
            
            if sum(sep_how == 'n')>0 %NMF    
                if initval
                    % background estimation:
                    %[out, bg(mm), bg_im]=backgroundoffset(dpixc);
                    [out, bg(mm), bg_im]=backgroundoffset(dpixc, 'no', 5, 20, 8); %empirical values...                  
                    
                    dvec_bg = ones(1, p.nx*p.ny);
                    %                     dvec_bg = p.offset*ones(1, p.nx*p.ny); %changed for
                    %                     offset 10...
                    
                    dvec_ind = squeeze(reshape(double(array2im(dpixc_ind)), p.nx*p.ny, 1, 2)); % vectors of resized images
%                     sum_dvec_ind = sum(dvec_ind, 1);
%                     dvec_ind = dvec_ind./repmat(sum_dvec_ind, p.nx*p.ny,1); %normlaized
                    
   
                    %                     hinit = [f*dvec_ind'; dvec_bg];       %original 'true' points + background
                    hinittmp = [dvec_ind'; dvec_bg];
                    sumh = sum(hinittmp,2);
                    hinit = hinittmp./repmat(sumh, 1, p.nx*p.ny); %normalized to 1
                    f = mean(dveccr(:)-bg(mm))/mean(mean(hinit(1:2,:))); %ration of the data/psf
                    %                     hinit = [rand(ncomp, p.nx*p.ny); dvec_bg];
                    blinkmatrand = rand(p.Nt, ncomp-1); %uniform random;
                    winit = [f*blinkmatrand, bg(mm)*sumh(ncomp)*ones(p.Nt,1)];             %random weights will be assigned to firts two and bg fixed

                    %                    [w{mm},h{mm}, wtrace{mm},htrace{mm}]=nmf_test(double(dveccr'),ncomp+1,1,winit,hinit, [3], [3]);
                    %                     [w{mm},h{mm}, wtrace,htrace,ddiv{mm}]=nmf_testconvD(double(dveccr'),ncomp+1,1,winit,hinit, [3], [3]);
                    [w{mm},h{mm}, wtrace,htrace,ddiv]=nmf_testconvD_normalH(double(dveccr'),ncomp,1,winit,hinit, [3], [3]);
                    qqq=[];
                else
                    [w{mm},h{mm}]=nmf(double(dveccr'),ncomp,1);
                end
                icapixNMF{mm} = reshape(h{mm}',p.nx,p.ny,ncomp);
            end
            
        end
        
        p.path_data = path_data;
        p.path_res = path_res;
        
        if savethis == 1
            fprintf('saving data \n');
            if ~(strcmp(path_res, path_data)) %not identical
                mkdir ([path_res namedir]);
                cd ([path_res namedir]);
            end
            %             save (p.namedir)
            save ([namedir '_separ'])
            writedata([],[],p,[namedir '_param'])
        end
    end
end

fprintf('\n')
