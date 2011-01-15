addpath('/home/liisa/080812/');
addpath('/home/rheintz/matlab/');

global para;

s=19;

% beadsAt=[180 504; 207 489; 253 236; 399 276; 262 54; 364 49; 393 168; 410 156; 300 139; 190 138; 376 282; 182 459; 333 449; 367 467; 415 453; 204 550; 175 579];

myres=para.res_resIm
fprintf('Please select bead coordinates by left clicking. Finish with right click\n');
beadsAt=dipgetcoords(100);
fprintf('%d beads selected\n',size(beadsAt,1));
beadsAt(end,:)=[];

[sumBdsHighRes,mybdsHRes,sumBdsSum,mybdsSum]=MeanFromCoord(myres,beadsAt,s,para.res_sumIm)

%myresSum=para.res_sumIm;
%[sumBdsSum,mybdsSum]=MeanFromCoord(myresSum,beadsAt,s)

%im = readtimeseries('beads1_000');
%myres=squeeze(sum(im(:,:,3:5)));
%[sumBdsRaw,mybdsRaw]=MeanFromCoord(myres,beadsAt,s)

% [params,res,fitted]=FitDataNDFast([10 10 0; 1000 s/2 s/2],sumBdsSum,300,'idiv')

[fittedSum,paramsSum,idiv,myfunct] = FitDataND('c(1)*exp(-((x{1}-c(2)).^2+(x{2}-c(3)).^2)/c(4))+c(5)',[max(sumBdsSum) 0 0 20 min(sumBdsSum)],sumBdsSum,1000);
sqrt(paramsSum(4))*para.res_resampledPixelsize(1) * 2*sqrt(log(2))  % to get FWHM
sumBdsSum-fittedSum

[fittedHRes,paramsHRes,idiv,myfunct] = FitDataND('c(1)*exp(-((x{1}-c(2)).^2+(x{2}-c(3)).^2)/c(4))+c(5)',[max(sumBdsHighRes) 0 0 20 min(sumBdsHighRes)],sumBdsHighRes,1000);
sqrt(paramsHRes(4))*para.res_resampledPixelsize(1)* 2*sqrt(log(2))  % to get FWHM
sumBdsHighRes-fittedHRes


