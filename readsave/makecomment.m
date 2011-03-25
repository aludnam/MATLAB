function makecomment(filename)
% makecomment(filename)
% Forces user to make commetn and stores it into the FILENAME
fid = fopen(filename,'a');
comment = input('Make some comments (line string):\n','s');
fprintf(fid, '\nComments:\n%s',comment);
fclose(fid);


