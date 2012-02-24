function makecomment(filename)
% makecomment(filename)
% Forces user to make commetn and stores it into the FILENAME
% If the file already exist, it will append the comment to the end of the
% file...
comment = input('Make some comments (line string):\n','s');
fid = fopen(filename,'a');
fprintf(fid, '\nComments (%s):\n%s',datestr(now),comment);
fclose(fid);