function printReadme(range)
% Prints readme files from the directories ['S' range(1)] unti ['S'
% range(end)]
% printReadme(range)
% Example: printReadme([320:325,333])


dirDefault = '~/project/data/qdots';
for ii=1:length(range)
    dirNow = ['S' num2str(range(ii))];
    if exist([dirDefault '/' dirNow '/readme'],'file')
        fprintf('Directory: %s\n',[dirDefault '/' dirNow]);
        type ([dirDefault '/' dirNow '/readme'])
        fprintf('===========================\n');
    end
end
        
    
    