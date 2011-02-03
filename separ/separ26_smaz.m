% evaluate simulated data - ICA and NMF
function separ(sep0, offset0, path_data, path_res, prename, niter, savethis, hint, sep_how)
% sep_how: string i - ica, n - nmf -> method for separatdion
if ~exist('hint', 'var')
    hint = 0;
end

if ~exist('sep_how', 'var')
    sep_how = 'in';
end

fprintf('Separation of components... \n')
for rr = 1: length(offset0)
    fprintf('\n%g: ',rr)
    for ll=1 : length(sep0)
        fprintf('.')
        p.namedir = [prename num2str(100*sep0(ll)) 'offset_' num2str(offset0(rr))];
        cd ([path_data p.namedir])
        for mm=1:niter
            
            namefile = [p.namedir '-iter_' num2str(mm)];
            load ([namefile '.mat'])
            %         ims(psf);
            %         SaveImageFULL('psf', 'pf');
            
            if sum(sep_how == 'i')>0 %ICA
                [icasig{mm}, A{mm}, W{mm}] = fastica (dveccr, 'numOfIC', 2, 'g', 'tanh');
                icapixICA{mm} = reshape(A{mm},32, 32, 2);
            end
            if sum(sep_how == 'n')>0 %NMF
                ncomp = 3; %number of components to be separated
                if hint
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
                    w{mm}=winit; h{mm}=hinit; %just testing...

                    %                    [w{mm},h{mm}, wtrace{mm},htrace{mm}]=nmf_test(double(dveccr'),ncomp+1,1,winit,hinit, [3], [3]);
                    %                     [w{mm},h{mm}, wtrace,htrace,ddiv{mm}]=nmf_testconvD(double(dveccr'),ncomp+1,1,winit,hinit, [3], [3]);
%                     [w{mm},h{mm}, wtrace,htrace,ddiv]=nmf_testconvD_normalH(double(dveccr'),ncomp,1,winit,hinit, [3], [3]);
                    
                else
                    [w{mm},h{mm}]=nmf(double(dveccr'),ncomp,1);
                end
                icapixNMF{mm} = reshape(h{mm}',32,32,ncomp);
            end
            %             imstiled(icapixICA{mm});
            %             SaveImageFULL([p.namedir 'ICA_' num2str(mm)], 'p');
            %             imstiled(icapixNMF{mm});
            %             SaveImageFULL([p.namedir 'NMF_' num2str(mm)], 'p');
            close all
            
        end
        
        p.path_data = path_data;
        p.path = path_res;
        
        if savethis == 1
            fprintf('saving data \n');
            if ~(strcmp(p.path, p.path_data)) %not identical
                mkdir ([p.path p.namedir]);
                cd ([p.path p.namedir]);
            end
            %             save (p.namedir)
            save ([p.namedir '_separ'])
            writedata([],[],p,[p.namedir '_param'])
        end
    end
end

fprintf('\n')
