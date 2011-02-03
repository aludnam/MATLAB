% evaluate simulated data - ICA and NMF
function separNMF(sep0, offset0, path_data, path_res, prename, niter, savethis, initval, sep_how)


if ~exist('initval', 'var')
    initval = 0;
end

if ~exist('sep_how', 'var')
    sep_how = 'in';
end
ncomp = 2; %number of components to be separated

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
                    blinkmatrand = rand(p.Nt, ncomp);
                    winit = [blinkmatrand,ones(p.Nt,1)];             %random weights will be assigned to firts two and bg fixed
                    
                    %                     [out, bg(mm), bg_im]=backgroundoffset(dpixc);
                    [out, bg(mm), bg_im]=backgroundoffset(dpixc, 'no', 5, 20, 8); %empirical values...
                    dvec_bg = bg(mm)*ones(1, p.nx*p.ny);
                    %                     dvec_bg = p.offset*ones(1, p.nx*p.ny); %changed for
                    %                     offset 10...
                    
                    dvec_ind = squeeze(reshape(double(array2im(dpixc_ind)), p.nx*p.ny, 1, 2)); % vectors of resized images
                    f = mean(dveccr(:))/mean(dvec_ind(:));
                    hinit = [f*dvec_ind'; dvec_bg];       %original 'true' points + background
                    %                     hinit = [rand(ncomp, p.nx*p.ny); dvec_bg];
                    ncomp = ncomp+1; %background added
                    %                    [w{mm},h{mm}, wtrace{mm},htrace{mm}]=nmf_test(double(dveccr'),ncomp+1,1,winit,hinit, [3], [3]);
                    %                     [w{mm},h{mm}, wtrace,htrace,ddiv{mm}]=nmf_testconvD(double(dveccr'),ncomp+1,1,winit,hinit, [3], [3]);
                    [w{mm},h{mm}]=nmf_testconvD(double(dveccr'),ncompbg,1,winit,hinit, [3], [3]);
                    
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
