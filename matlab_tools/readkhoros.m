%READKHOROS   Read image in khoros file format from disk
%
% SYNOPSIS:
%  readkhoros(filename)
%
% PARAMETERS:
%  filename:    string with name of file.

% (C) Copyright 2004           Department of Molecular Biology
%     All rights reserved      Max-Planck-Institute for Biophysical Chemistry
%                              Am Fassberg 11, 37077 G"ottingen
%                              Germany
%
% Bernd Rieger & Rainer Heintzmann, June, 2004.


function out=readkhoros(fn)

fid = fopen(fn);
if (fid < 0) 
  error(['ReadKhoros: File ' fn ' not found on disc.']);
end
[dims, dtype] = readkhoros_info(fn);
if isempty(dtype)
    fprintf('WARNING: khoros file unreadable. Trying conversion to new format\n');
    command=sprintf('/usr/local/KhorosInstall/bin/kcptoval -val -i %s -o %s',fn,fn)
    system(command);
    [dims, dtype] = readkhoros_info(fn);
end    

status = fseek(fid,dims(6),-1);
if status <0
   error(ferror(fid));
end
dims(end) = [];
ii=5; % go backwards in dimensions
while dims(ii) ==1 & ii > 1 %khoros images max 5D
   dims(ii)=[]; %remove trailing ones
   ii=ii-1;     %but keep in the middles ones
end

if findstr(dtype,'scomplex')
    [out,count]=fread(fid,[prod(dims)]*2,'single => single');
%error('Sorry no complex value reading.');
elseif findstr(dtype,'dcomplex')
    [out,count]=fread(fid,[prod(dims)]*2,'double');
else
    [out,count]=fread(fid,[prod(dims)],[dtype '=>' dtype]);
end

if count==0
   error('Could not read file.');
end

if findstr(dtype,'complex')
    out=reshape(out,[2 size(out,1)/2]);
    out=out(1,:)+i*out(2,:);
end

out = dip_image(out,dtype);
%out = squeeze(reshape(out,[dims']));
out = reshape(out,[dims']);

fclose(fid);
