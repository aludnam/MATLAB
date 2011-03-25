function copyscript(source, destination, source_app, destination_app)
% copyscript(source, destination, source_app, destination_app)
% name format [source '_script' source_app '.m']
if ~exist('source_app','var')
    source_app = [];
end
if ~exist('destination_app','var')
    destination_app = [];
end

dest_dir_pref='~/project/data/qdots/';
dest_dir = [dest_dir_pref destination];
mkdirsafe(dest_dir);

fprintf('Making directory: %s\n',dest_dir)
app = '_script';
source_all = [source app source_app '.m'];
dest_all = [dest_dir '/' destination app destination_app '.m'];
copyfile (source_all, dest_all)
fprintf('Copying file %s to %s\n' , source_all, dest_all);
if exist('params.m', 'file')
    copyfile ('params.m', [dest_dir '/params.m'])
    fprintf('Copying parameters\n');
end

cd (dest_dir)
makecomment('readme');