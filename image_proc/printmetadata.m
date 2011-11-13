function printmetadata(pathname, metadataFileName)
% printmetadata(pathname)
% Extracts Time, Comment and Exposure from metadataFileName. (default: metadata.txt)
if ~exist('metadataFileName','var')
    metadataFileName = 'metadata.txt';
end
fprintf('Directory %s\t (Figure %2.0f):\n',pathname,gcf);    
system (['cat ' pathname '/' metadataFileName '| grep \"Time\": -m 1']);   
system (['cat ' pathname '/' metadataFileName '| grep Comment']);    
system (['cat ' pathname '/' metadataFileName '| grep Exposure -m 1']);    