function mfprintf(fid, varargin)
% mfprintf(fid, varargin)
% As fprintf but allows fid to be a vector. 
% For example:
% fid = fopen('foo.txt','a');
% mfprintf([fid 1], 'Hello %s\m','world')
% fclose(fid);
% prints 'Hello world' into file foo.txt and on the screen at the same
% time.

arrayfun(@(fid) fprintf(fid, varargin{:}), fid,'UniformOutput',0);
