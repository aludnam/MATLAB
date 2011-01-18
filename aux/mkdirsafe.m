function [success,message,messageid] = mkdirsafe(name_directory)
% [success,message,messageid] = mkdirsafe(name_directory)
% Makes direcotry safely. if it already exist gives error

[success,message,messageid] = mkdir(name_directory);

if ~isempty(message)
    error(message)
end
