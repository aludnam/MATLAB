function writedata (x, y, p, filename, commentStr, dataonly)
% WRITEDATA creates txt file with parameters and computed values
% writedata (x, y, p, filename, commentStr)
% x - values on x-axis
% y - values on y-axis
% p - parameters used for computation
% filename - name of the txt file
% commentStr - comment string
% dataonly - binary. if dataonly == 1 writes only numeric ASCI of data x
% and y (doesn't write p, or comment string), dataonly==0 (default) - adds comment
% string

fprintf('Writing into file: %s\n', filename);
fid = fopen([filename '.dat'],'wt');
if nargin<6
    dataonly = 0;
end

if dataonly==0
    
    if nargin>4
        fprintf(fid,'%s \n', commentStr);
    end
    
    fprintf(fid,'Computation: %s \n\n', datestr(now));
    
    fprintf(fid,'Parameters:\n');
    
    if ~isempty(p)
        if isstruct(p)
            fn =fieldnames(p);
            for ii = 1: length(fn)
                name = fn{ii};
                value = getfield(p,fn{ii});
                if ischar (value)
                    fprintf (fid, [name ' = ''' value '''\n']);
                else 
                    if isempty(value)
                        fprintf (fid, [name ' = [ ]\n']);                        
                    elseif length(value)>1
                        fprintf (fid, [name ' = [']);
                        fprintf (fid, '%g ', value);
                        fprintf (fid, ']\n');                        
                    else
                        fprintf (fid, [name ' = %g \n'], value);
                    end
                end
            end
        else
            for ii = 1: length(p)
                fprintf (fid, ['p' num2str(ii) ' = %g \n'], p(ii));
            end
            
        end
    end    
end

if or(~isempty(x), ~isempty(y))
    fprintf (fid,'\nX \t\t\t Y \n');
    fprintf (fid,'------------------\n');
    for ii=1:length(y)
        fprintf (fid,'%-8.3f \t %-8.3f \n', x(ii), y(ii));
    end
end
fclose(fid);
