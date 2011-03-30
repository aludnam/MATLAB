function [logfile, fid] = initlogfile
% [logfile, fid] = initlogfile
% Initializes log fiel in the current directory and gives FID (vector) for
% mfprintf fucntion.
logfile = 'evaluation.log';
fidlog = fopen(logfile,'wt');
fprintf(fidlog,'Computation: %s\n\n', datestr(now));
fid = [1 fidlog];
