function copyscript(source, destination, source_app, destination_app)
% copyscript(source, destination, source_app, destination_app)
% name format [source '_script' source_app '.m']
% Examples:
% copyscript('S333', 'S334') will make a folder S334 and copy the script S333_script.m and the
% params.m from the folder S333 with correct name (S334/S334_script.m).
% copyscript([],'S335') will make a folder S335 and copy hte Default script
% (Default/Default_script.m) and parameters (Default/params.m) to this
% folder with correct name (S335/S335_script.m).

if ~exist('source_app','var')
    source_app = [];
end
if ~exist('destination_app','var')
    destination_app = [];
end

dest_dir_pref='~/project/data/qdots';
dest_dir = [dest_dir_pref '/' destination];
fprintf('Making directory: %s\n',dest_dir)
mkdirsafe(dest_dir);

if isempty(source)
    source_dir = [dest_dir_pref '/Default'];
    source = 'Default';
    fprintf('Default directory will be used as the source.\n')
else 
    source_dir = cd; 
end

app = '_script';
source_all = [source_dir '/' source app source_app '.m'];
dest_all = [dest_dir '/' destination app destination_app '.m'];
copyfile (source_all, dest_all)
fprintf('Copying file %s to %s\n' , source_all, dest_all);
if exist([source_dir '/params.m'],'file')
    copyfile ([source_dir '/params.m'], [dest_dir '/params.m'])
    fprintf('Copying parameters\n');
end

cd (dest_dir)
makecomment('readme');