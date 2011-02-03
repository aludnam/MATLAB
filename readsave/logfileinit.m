function fid = logfileinit(logfile)
% fid = logfileinit(logfile)
fid = fopen(logfile,'wt');
fprintf(fid,'Computation: %s\n\n', datestr(now));

