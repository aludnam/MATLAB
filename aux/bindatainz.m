function data_out = bindatainz(data)
% data_out = bindatainz(data)
sz= size(data, 3);
ixo=[1:2:sz];
ixe=[2:2:sz];
data_out = data(:,:,ixe)+data(:,:,ixo(1:length(ixe)));
