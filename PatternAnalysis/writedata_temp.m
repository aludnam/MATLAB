function writedata (x, y, yMax, yMin, p, filename, commentStr)
% WRITEDATA creates txt file with parameters and computed values
% writedata (x, y, yMax, yMin, p, filename, commentStr)
% x - values on x-axis
% y - values on y-axis
% p - parameters used for computation 
% filename - name of the txt file
% commentStr - comment string
fid = fopen([filename '.txt'],'wt');

if nargin>6
    fprintf(fid,'%s \n', commentStr);
end
    
fprintf(fid,'Computation: %s \n\n', datestr(now));

fprintf(fid,'Parameters:\n');


fn =fieldnames(p);
for ii = 1: length(fn)
    name = fn{ii};
    value = getfield(p,fn{ii});
    if ischar (value)
        fprintf (fid, [name ' = ' value '\n']);
    else
        fprintf (fid, [name ' = %g \n'], value);
    end   
end

fprintf (fid,'\nX \t\t\t Y \n');
fprintf (fid,'------------------\n');

for ii=1:length(y)
    fprintf (fid,'%-8.3f \t %-8.3f \n', x(ii), y(ii));
end

fclose(fid);


%envelopes 

if ~isempty(yMax)
    fid = fopen([filename '_envelope_Max.txt'],'wt');
    if nargin>6
        fprintf(fid,'%s \n', commentStr);
    end
    fprintf(fid,'Computation: %s \n\n', datestr(now));
    fprintf(fid,'Maximum of %g simulations of Poisson process\n\n', p.Nsimul);
    fprintf (fid,'x \t\t\t max y \n');
    fprintf (fid,'------------------\n');    
    for ii=1:length(y)
        fprintf (fid,'%-8.3f \t %-8.3f \n', x(ii), yMax(ii));
    end
    fclose(fid);
end


if ~isempty(yMin)
    fid = fopen([filename '_envelope_Min.txt'],'wt');
    if nargin>6
        fprintf(fid,'%s \n', commentStr);
    end
    fprintf(fid,'Computation: %s \n\n', datestr(now));
    fprintf(fid,'Minimum of %g simulations of Poisson process\n\n', p.Nsimul);
    fprintf (fid,'xK \t\t\t min K \n');
    fprintf (fid,'------------------\n');    
    for ii=1:length(y)
        fprintf (fid,'%-8.3f \t %-8.3f \n', x(ii), yMin(ii));
    end
    fclose(fid);
end

