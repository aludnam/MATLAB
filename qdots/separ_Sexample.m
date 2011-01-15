% evaluate simulated data - NMF
% V ~ WH
% V -> N_pix x N_t
% W -> N_pix x N_comp - xth pixel of the ith components
% H -> N_copm x N_t - contribution of the i-th component in the time t

function separ_Sexample(sep0, offset0, path_data, path_res, prename, niter, avg_num, savethis)

ncomp = 3; %number of components to be separated + background
randh=rand(2,1,4); %random coeeficients for hinit

for rr = 1: length(offset0)
    for ll=1 : length(sep0)
        for qq=1 : length(avg_num)
            namedir = [prename num2str(100*sep0(ll)) 'offset_' num2str(offset0(rr))];
            cd ([path_data namedir])
            for mm=1:niter
                %reads the first...
                namefile = [namedir '-iter_' num2str(mm)];
                load ([namefile '.mat'])
                %cat the folowing,,,
                p.catitervec=[1];
                [dpixc, dveccr, blinkmat, p] = catsimul(namedir, p.catitervec);
                
                [out, bg(mm), bg_im]=backgroundoffset(dpixc, 'no', 5, 20, 8); %empirical values...               
                dvec_ind = squeeze(reshape(double(array2im(dpixc_ind)), p.nx*p.ny, 1, 2)); % vectors of resized images                                
                dvec_bg = ones(p.nx*p.ny, 1);
                winittmp = [dvec_ind, dvec_bg];
                sumw = sum(winittmp,1);
                winit = winittmp./repmat(sumw, p.nx*p.ny, 1); %normalized to 1                
                p.Nt=1;
                p.meanblinkmat=mean(mean(blinkmat,2));
                p.htrue=[0.3; 0.7];
                htruef=p.meanblinkmat*p.htrue;
                
                for kk=1:4
                    hinit = [p.meanblinkmat*randh(:,:,kk); p.offset*sumw(ncomp)*ones(1, p.Nt)];
                    %hinit = [htruef; p.offset*sumw(ncomp)*ones(1, p.Nt)];
                    dveccr=dvec_ind*htruef+p.offset;
                    dpixcd = dip_image(reshape(dveccr,p.nx,p.ny));
                    dpixcdn_tmp = newim(dpixcd);
                    for oo=1:avg_num(qq)
                        dpixcdn_tmp = dpixcdn_tmp + noise(dpixcd,'poisson');
                    end
                    if avg_num(qq)==0
                        dpixcdn_tmp=dpixcd; %noise free
                    end
                    dpixcdn=dpixcdn_tmp/max(avg_num(qq),1); %average of (avg_num(qq))X realization...
                    dveccr=reshape(double(dpixcdn),p.nx*p.ny,1);
                    %[c,w{mm},h{mm}, X1,X2, dhr, minXr, miXvalr, mhdr, htrace, p]
                    [c,w{mm},h{mm}, X1,X2, dh, minX, miXval, mhd, htr,htrace, p]=nmf_S51(double(dveccr),ncomp,1,winit,hinit, [3], [3],p);
                    res{kk}=struct('c',c,'w',w,'h',h,'X1',X1,'X2',X2,'dh',dh,'minX',minX,'miXval',miXval, 'mhd',mhd, 'htr',htr, 'htrace',htrace,'p',p);
                end
                
                
                %hinit = [f*blinkmat./repmat(mean(blinkmat,2),1,size(blinkmat,2)); bg(mm)*sumw(ncomp)*ones(1, p.Nt)];             %random weights will be assigned to firts two and bg fixed
                %[c,w{mm},h{mm}, X1,X2, dht, minXt, miXvalt, mhdt, htrace,p]=nmf_S41(double(dveccr),ncomp,1,winit,hinit, [3], [3],p);
                qqq=[];
                icapixNMF{mm} = reshape(w{mm},p.nx,p.ny,ncomp);
                
                
            end
            p.avg_num = avg_num(qq);
            p.path_data = path_data;
            p.path_res = path_res;
            
            if savethis == 1
                fprintf('saving data \n');
                if ~(strcmp(path_res, path_data)) %not identical
                    mkdir ([path_res namedir]);
                    cd ([path_res namedir]);
                end
                qqq=[];
                %             save (p.namedir)
                %save ([namedir '_separMulti'],save ([namedir
                %'_separMulti'],'p', 'X1','X2','dhr','minXr',
                %'miXvalr','dht'))
                save ([namedir '_DdviMap_avg' num2str(avg_num(qq))],'res')
                writedata([],[],p,[namedir '_avg' num2str(avg_num(qq)) '_param'])
            end
        end
    end
end

fprintf('\n')
