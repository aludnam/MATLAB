[afilename, apathname] = uigetfile('*.*', 'Input for Beads');
filename=fullfile(apathname,afilename);
%filename='E:\2009-05-12_GlassesBeads\tif\tryout\163m_glass40.tif'

s=40;
%a=ZEN_LSMLoad(filename);
a=readtimeseries(filename);

a=a{1};
%%Original version 20/05/2009
% ma=max(a{1},[],3);
% ExtractResults(ma);

%%Marie, 21/05/2009
% ma=max(a{1},[],3); 
[pxy,bestslice]=ExtractCoordinates(a,s);

number=size(bestslice,2);
for b=1:number
    slice=squeeze(extract(a,[s s 1],[pxy(b,:) bestslice(b)])); 

    FWHM(b)=FitPSF(slice); 
end
sumFWHM=0;
for b=1:number
   sumFWHM=sumFWHM+FWHM(b); 
end
fprintf('The average FWHM in pxl of the selected beads is\n');
averageFWHM=sumFWHM/number
fprintf('The average FWHM in nm of the selected beads is\n');
ConvertZENpxl_nm(averageFWHM)