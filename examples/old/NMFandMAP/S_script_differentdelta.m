cleardataquestion   % clearing the data
savethis = 1;       % saves the results
params              % parameters for the computitation

delta_vec = [.1 .3 .5 .8 1]; %all...
ncomp_vec = [2];
for nc_ix = 1:length(ncomp_vec)
    peval.ncomp = ncomp_vec(nc_ix); %number of components without the background component....
    for delta = delta_vec
        % data location for each computaution:
        peval.data_path = [getenv('HOME') '/project/data/qdots/D2ptBMrandOblique' '/'];
        peval.data_dir= ['D_N2_delta' num2str(delta*100) '_NoNoise'];
        peval.data_name=peval.data_dir;
        
        % results name:
        peval.res_nameappendix = ['_nc' num2str(peval.ncomp)];
        peval.res_dir = [peval.data_dir peval.res_nameappendix];
        peval.res_name = peval.res_dir;
        
        %reads data:
        [dpixc, dpixc_ind, blinkmat, peval, p] = readSimulSimple(peval);
        peval.bg = peval.offset; %cheating here
        
        for iteration=1:10 % several evaluations of the same data
            for nn = 1 : peval.nrestartH
                % initialization:
                winit_pix = init_w_meandata1res;
                hinit = rand(peval.ncomp, size(dpixc,3));
                %nmf separation:
                [res, peval] = separcomp(dpixc, peval, winit_pix, hinit);
                
                %saving results
                if savethis
                    saveresultscompleteSimple
                end
            end
            
            % MAP fiting with nmf initialization
            peval.filenamebase = ['iter' num2str(iteration)];
            % peval.sigmaPSF = 1.3; %sigma of the PSF gauss approx (in pixels)
            if savethis; cd (peval.res_dir); end
            [ihs, x_mu, y_mu, sig] = plotreswh2(res, peval, dpixc, p, savethis, 1, 0, 1, 1,0);
            if savethis; cd ..; end
            % initialization form nmf ressults:            
            x1in=reshape(res.h(1:peval.ncomp,:), 1,peval.ncomp*peval.nt);
            x2in=-1+[x_mu, y_mu];
            [x1in,x2in]=initx1x2(Hinit,Winit,x1,x2,x1rand,x2rand);
            Wres=x2in;
            Hres=x1in;
            pointlogWall=[]; flog=[]; pointlogHall=[];
            for kk=1:peval.nIterAlter %alternating optimization
                
                % W fixed, optimizing H (intensities)
                [Hres, optionsH, flogH, pointlogH] = conjgrad(peval.Hupdate, Hres, optionsH, ['grad' peval.Hupdate],reshape(dpixc, peval.nx*peval.ny, peval.nt), peval.sigmaPSF, peval.alpha, peval.beta, peval, Wres(1:peval.ncomp), Wres(peval.ncomp+1:end));
                flog=[flog; flogH];
                pointlogHall=[pointlogHall; pointlogH];
                % H fixed, optimizing W (positions)
                [Wres, optionsW, flogW, pointlogW] = conjgrad(peval.Wupdate, Wres, optionsW, ['grad' peval.Wupdate],reshape(dpixc, peval.nx*peval.ny, peval.nt), peval.sigmaPSF, peval.alpha, peval.beta, peval, Hres);
                flog=[flog; flogW];
                pointlogWall=[pointlogWall; pointlogW];
                                
                
                peval.Wres=Wres;
                %% plotting & saves results
                if savethis; cd (peval.res_dir); end
                showimage=0;
                res.dpixc=dpixc;
                plotloglikGaP;
                if savethis; cd ..; end
                
            end
        end
    end
end
