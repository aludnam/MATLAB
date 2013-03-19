% SaveImageFULL (nameImage, saveFormat, handleImage)
function SaveImageFULL (nameImage, saveFormat, handleImage)
if nargin<2
    saveFormat = 'efpd';
end
if nargin<3 
    handleImage = gcf;
end
fprintf('Saving image: %s (format: ''%s'')\n', nameImage, saveFormat);
pth=cd;
fprintf('(path: %s)\n',pth)
if sum(saveFormat == 'e')>0 saveas(handleImage, nameImage,'epsc'); end
if sum(saveFormat == 'p')>0 saveas(handleImage, nameImage,'png'); end
if sum(saveFormat == 'f')>0 saveas(handleImage, nameImage,'fig'); end
if sum(saveFormat == 'd')>0 saveas(handleImage, nameImage,'pdf'); end

