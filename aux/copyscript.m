function copyscript(source, source_app, destination, destination_app)
% copyscript(source, source_app, destination, destination_app)

dest_dir = ['../' destination];
[success,message,messageid] = mkdir(dest_dir);

if ~isempty(message)
    error(message)
end

fprintf('Making directory: %s\n',dest_dir)
app = '_script';
source_all = [source app source_app '.m'];
dest_all = [dest_dir '/' destination app destination_app '.m'];
copyfile (source_all, dest_all)
cd (dest_dir)
fprintf('Copying file %s to %s\n' , source_all, dest_all);
if exist('params.m', 'file')
    copyfile ('params.m', [dest_dir '/params.m'])
    fprintf('Copying parameters\n');
end

