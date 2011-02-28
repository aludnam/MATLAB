function [success,message,messageid] = mkdirsafe(name_directory)
% [success,message,messageid] = mkdirsafe(name_directory)
% Makes direcotry safely. Gives error if it already exist.

[success,message,messageid] = mkdir(name_directory);

if ~isempty(message)
    error(message)
end
