R = input ('Clear data? [y/n] \n','s');
if strcmp(R,'y')
    clear
end
savethis = 1;
peval.res_path = [cd '/'];

%simulation parameters:
peval.N = 5;
peval.offset = 100;
peval.niter_vec = [1];
peval.avg_num_vec = [1]; %number of averages os samples (reducing Poisson noise) avg_num=0 --> noise free image...
peval.catiter_vec=[1]; %cat iterations of simulations - to make data longer
peval.ddterm = eps;
peval.maxiter = 500;


%evaluation parameters
ncomp_vec = peval.N;
peval.w_fixvec = []; %figing only added backround component
peval.h_fixvec = [];

peval.maxh = 10;
peval.maxw = 10;
peval.method = 'nmf_classic';
peval.nrestartH = 5; %number of restarts of H
% ncomp_vec = [1 5 10 15 18 20 22 25 28 30 32 35 40];
%fitting part parameters
peval.Witer=3;          %number of updates of W
peval.Hiter=20;         %number of updates of H
peval.nIterAlter = 10;  %number of cycles of updates W and H

peval.Wupdate='loglikGaPExpHfix';
peval.Hupdate='loglikGaPExpWfix';

% peval.ROIx = [1,9];
% peval.ROIy = [1,9];
% peval.ROIz = [1,size(dpixc,3)];
%
% ROIx=peval.ROIx(1):peval.ROIx(2);
% ROIy=peval.ROIy(1):peval.ROIy(2);
% ROIz=peval.ROIz(1):peval.ROIz(2);
ROIx=[];
ROIy=[];
ROIz=[];




iter = peval.niter_vec(1);
avg = peval.avg_num_vec(1);
%%%%%%%%%%%%%%%%%%%%%%%%
% peval.filenamebase = ['ExpBlinkmatSigmaPSF' num2str(peval.sigmaPSF*100)
% '-ExpPrior'];
optionsW = zeros(1,18);
optionsW(1)=1;                  %to display error values
optionsW(7)=1;
optionsW(9)=0;                   %to check gradient
optionsW(14)=peval.Witer;         %maximum number of iterations

optionsH=optionsW;
optionsH(14)=peval.Hiter;

Hinitvec{1}='trueH'; Hinitvec{2}='wrongH'; Winitvec{1}='trueW'; Winitvec{2}='wrongW';
%fitting part parameters
%%%%%%%%%%%%%%%%%%%%%%%%
% delta_vec=[.6 1.0 1.2 1.4 1.6 1.8 2.0];
% delta_vec=[.6 .8 1.0 1.2 1.4 1.6 1.8 2.0];
% delta_vec=[.4 .6 .8 1.0 1.2 1.4 1.6 1.8 2.0];
delta_vec=[1.2 1.4 1.6 1.8 2.0];
for delta = delta_vec
    for iteration=2:10
        peval.data_dir= ['D_N5_delta' num2str(delta*100)];
        peval.data_name=peval.data_dir;
        peval.data_path = [getenv('HOME') '/project/data/qdots/D5ptBMrandStrong' '/'];
        
        
        for nc_ix = 1:length(ncomp_vec)
            peval.ncomp = ncomp_vec(nc_ix); %number of components without the background component....
            peval.addnamefileres = ['_nc' num2str(peval.ncomp)];
            peval.res_dir = [peval.data_dir peval.addnamefileres];
            peval.res_name = peval.res_dir;
            [dpixc, dpixc_ind, blinkmat, peval, p] = readSimulSimple(peval);
            if ~isempty(ROIx*ROIy*ROIz) %only when all are not empty
                roicutdata %cut the data to the ROI
            end
            
            for nn = 1 : peval.nrestartH
                if nn==1
                    a=mean(dpixc,3);
                    an=a/sum(a(:));
                    winit_pix(:,:,1:peval.ncomp) = repmat(an,[1,1,peval.ncomp]);
                else
                    ca
                    winit_pix(:,:,1:peval.ncomp) = reshape(res.w(:,1:peval.ncomp), peval.nx, peval.ny, peval.ncomp);
                end
                %                 hinit = [];
                hinit = rand(peval.ncomp, size(dpixc,3));
                peval.bg = peval.offset;
                tstart = tic;
                [res, peval] = separcomp(dpixc, peval, winit_pix, hinit);
                peval.elapsedtime = toc(tstart);

                %saving results
                if savethis>0
                    saveresultscompleteSimple
                end
                
            end
            
            %form NMF:
            peval.sigmaPSF = p.sigmapix; %sigma of the PSF gauss approx (in pixels)
            peval.alpha = 1;        %exponentional prior on H
            peval.beta=1000;
%             peval.filenamebase = ['FitSigmaPSF' num2str(round(peval.sigmaPSF*100)) '-' peval.Wupdate];
            peval.filenamebase = ['iter' num2str(iteration)];
            % peval.sigmaPSF = 1.3; %sigma of the PSF gauss approx (in pixels)
            if savethis; cd (peval.res_dir); end
            [ihs, x_mu, y_mu, sig] = plotreswh2(res, peval, dpixc, p, savethis, 1, 0, 1, 1,0);
            if savethis; cd ..; end
%             % resnmf=load('../S170/D5ptBMexp750PSFgauss_N5_offset100_iter1_avg1_nc5/D5ptBMexp750PSFgauss_N5_offset100_iter1_avg1_nc5.mat');
            %x1=reshape(blinkmat, 1,peval.ncomp*peval.nt);
            x1=[];
            x2=[p.x_vec, p.y_vec];
            x1rand=reshape(res.h(1:peval.ncomp,:), 1,peval.ncomp*peval.nt);
            x2rand=-1+[x_mu, y_mu];
            % x2rand=x2+5*(rand(1, 4))-2.5;
            
            for indH=2
                for indW=2
                    %initialization:
                    Hinit=Hinitvec{indH};
                    Winit=Winitvec{indW};
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
                        margin=.2;
%                         figure(2)
%                         hold on
%                         scatter(p.x_vec, p.y_vec)
%                         set(gca,'YDir','reverse')
%                         for ii=1:peval.ncomp
%                             scatter(Wres(:,end-peval.ncomp-ii+1), Wres(:,end-ii+1),'x')
%                         end
%                         % scatter(Wres(end-2), Wres(end),'dr')
%                         % scatter(Wres(end-3), Wres(end-1),'xr')
%                         grid on
%                         l=max(p.x_vec)-min(p.x_vec)+2*margin;
%                         
%                       
%                         xlim([min(p.x_vec)-margin max(p.x_vec)+margin])
%                         ylim([mean(p.y_vec)-l/2 mean(p.y_vec)+l/2])
                        
                    end
                    peval.Wres=Wres;
                    %% plotting & saves results
                    if savethis; cd (peval.res_dir); end
                    showimage=0;                   
                    plotloglikGaP;
                    if savethis; cd ..; end
                end
            end
        end
    end
end
% [Wxk,Hkt,centers,Wxkpix]=reshapeGaP(res.hvec,res.cxcy,peval);
% Vxtpixbg=reshape(Wxk*Hkt,peval.nx,peval.ny,peval.nt)+peval.bg;
% resid=(Vxtpixbg-dpixc);
% resid_norm=resid./Vxtpixbg;
% % [Z,H,T,perm] = dendrogramfromtimecorrelation(resid_norm,'complete'); %'average');
% [Z,H,T,perm] = dendrogram_subtreepixels(resid_norm,'complete', p , x_mu, y_mu, savethis)
