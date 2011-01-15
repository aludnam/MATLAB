% [normalized,TotalBleach,params]=EvalFRAP('filename',whichroi,maxframe,prebleach): Evaluates LSM FRAP data by requesting various ROIs and normalising and plotting and fitting the results.
% filename : LSM file to process
% whichroi (optional) : Name the ROI number (if multiples are present) where the bleach was.
% maxframe : limit the analysis to a maximal number of frames (default = -1 = all frames)
% prebleach : number of prebleach images (default = 5)
% normalized : normalised FRAP data
% TotalBleach : Curve of total decay during experiment
% params : Fit results [BleachDepth, Recovery, t2]
% Written for Fernanda
%
function [normalized,timepoints,TotalBleach,params]=EvalFrap(filename,whichroi,maxframe,prebleach)
if nargin < 2
    whichroi=-1; % all rois are valid
end
if nargin < 3
    maxframe=-1; % all rois are valid
end
if nargin < 4
    prebleach=5; % number of pre bleach images
end
[data,ts,myroi,myroi2]=kLoadLSM(filename);
if maxframe > 0
    data=data(:,:,:,0:maxframe-1);
    ts=ts(0:maxframe-1);
end
% 'This is the bleach ROI'
roi=cat(3,myroi,myroi2);
if whichroi > 0
    roi=roi(:,:,whichroi);
end
hasroi=sum(roi,[],[1 2]) > 0;
if sum(hasroi) ~= 1
    roi
    error('Data not exactly one ROI');
end
roi=sum(roi,[],3) > 0;

img = mean(data,[],4)
fprintf('Select a ROI to evaluate\n');
[evalROI,v] =diproi(gcf);
evalROI=repmat(evalROI,[1 1 size(data,4)]);

dat2=croptomask(squeeze(data),evalROI)
sz=size(data,4);
clear data

roi2c=croptomask(roi,squeeze(evalROI(:,:,0)));
clear roi
BleachROI=repmat(roi2c,[1 1 sz]);
% Align
fprintf('Aligning Data\n');
dat2c=kcorrelator(dat2,{'p',1;'fixplane',prebleach-1});  % Use plane 4 as reference plane

mImg2 = squeeze(mean(dat2c,[],3));
[datBg, bg, bg_im]=backgroundoffset(mImg2);
bg_im
fprintf('Background is %d\nSelect Reference ROI\n',bg)
mImg2
RefROI=diproi(gcf);
RefROI=repmat(RefROI,[1 1 sz]);

dat2cb=dat2c-bg;

bleach=squeeze(mean(dat2cb,BleachROI,[1 2]));
preBl = mean(bleach(0:prebleach-1));
ref=squeeze(mean(dat2cb,RefROI,[1 2]));
preRef = mean(ref(0:prebleach-1));

TotalBleach=ref/preRef

normalized = (bleach/preBl) * preRef /ref;
tsn=ts-ts(prebleach);
timepoints=tsn;  % for return
plot(tsn,normalized)
hold on

recNor=squeeze(normalized(prebleach:end));
tNor=squeeze(tsn(prebleach:end));

tFitStart=tNor(end)/4;

[fitted,params,myfunct,msevalue] = FitNlinData('c(1)+(1-c(1))*c(2)*x./(c(3)+x)',double([recNor(1) recNor(end)-recNor(1) tFitStart]),recNor,tNor);

plot(tNor,fitted,'red')
legend(['Data Bacgkround ' num2str(bg)],['BleachDepth: ' num2str(params(1)) ', Recovery: ' num2str(params(2)) ', t2: ' num2str(params(3)) ' s'])
xlabel('t [s]')
ylabel('normalized recovery')

%dsjkshgdk

% save all that needs to be saved
filebase=filename(1:end-4);
print([filebase '_FRAP.eps'],'-depsc')
print([filebase '_FRAP.jpg'],'-djpeg90')

save([filebase '_fit.txt'],'params','-ASCII')
tosave=double([tNor;fitted])';
save([filebase '_fit_curve.txt'],'tosave','-ASCII')
tosave=double([tsn;normalized])';
save([filebase '_curve.txt'],'tosave','-ASCII')


tiffwrite([filebase '_BgROI.tif'],bg_im,'yes')
tiffwrite([filebase '_BgIm.tif'],mImg2,'yes')

