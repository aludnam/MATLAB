function imgs = readspefile(inputfile)
% imgs = readspefile(inputfile)
% Reads SPE file. This is copied from CSSTORM.m

fid = fopen(inputfile, 'r', 'l');
header = fread(fid,2050,'*uint16');

% parsing .spe header
x_dim = double(header(22));
y_dim = double(header(329));
z_dim = double(header(724));
img_mode = header(55);
if img_mode == 0
    imgs = fread(fid, inf, 'float32');
elseif img_mode == 1
    imgs = fread(fid, inf, 'uint32');
elseif img_mode == 2
    imgs = fread(fid, inf, 'int16');
elseif img_mode == 3
    imgs = fread(fid, inf, 'uint16');
end

fclose(fid);

imgs = reshape(imgs, [x_dim y_dim z_dim]);