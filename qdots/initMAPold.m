function [x,y, hinit] = initMAP(dpixc, peval)

peval.meandata = mean(dpixc(:));
image=mean(dpixc,3);         
[x,y]=sample2d(image,peval.ncomp);

hinit = rand(peval.ncomp, size(dpixc,3));
winit_pix_tmp=gauss2dmultislice([peval.nx, peval.ny, peval.ncomp], [x',y']+1, peval.sigmaPSF, 1);
winit_pix = normalizeslices(winit_pix_tmp);
[winit, hinit] = initwh(winit_pix, hinit, peval); %initialization of w and h