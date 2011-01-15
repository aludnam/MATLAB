%WRITEKHOROS   Write image in khoros file format
%
% SYNOPSIS:
%  writekhoros(image_in,filename)
%
% PARAMETERS:
%  image_in:    image to write to file.
%  filename:    string with name of file, optionally with path and extension.

% (C) Copyright 2004           Department of Molecular Biology
%     All rights reserved      Max-Planck-Institute for Biophysical Chemistry
%                              Am Fassberg 11, 37077 G"ottingen
%                              Germany
%
% Bernd Rieger & Rainer Heintzmann, June, 2004.


function writekhoros(in,fn)


fid = fopen(fn,'a');
if fid <0
   error(['Could not open file ' fn  ' to append data.']);
end
fclose(fid);
   
if isa(in,'dip_image')
   dim = [1 1 1 1 1];
   sz=size(in);
   if length(sz)>5
      error('array max 5d.');
   end
   dim(1:length(sz))=sz;
   writekhoros_info(fn,dim,datatype(in));
   dipio_appendrawdata(in,fn);
else
   dim = [1 1 1 1 1];
   sz=size(in);
   if length(sz)>5
      error('array max 5d.');
   end
   dim(1:length(sz))=sz;
   writekhoros_info(fn,dim,'double');
   fid = fopen(fn,'a');
   fwrite(fid,in,'double');
   fclose(fid);
end

