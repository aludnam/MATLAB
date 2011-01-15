function  plotHtiled_all(pathres, sep, offset)
% cd ~project/data/qdots/S50
% sep=[10 20 50];
% offset=[100 1000];
cd (pathres)
sufx = {'avg1', 'avg2', 'avg10', 'avg50', 'avg0'};
% col = 'rbgkm';
col = get (gca,'colororder');
for is=1:size(sep,2)
    for io=1:size(offset,2)
        figure
        for ii=1:size(sufx,2)
            namebase=['S44_sep_' num2str(sep(is)) 'offset_' num2str(offset(io))];
            nameload=[namebase '/' namebase '_DdviMap_' sufx{ii} '.mat'];
            load (nameload)
            plotHtiled(res,1, col(ii,:))
%             SaveImageFULL([namebase '/h_iter_' sufx{ii}], 'epf')
        end
        SaveImageFULL([namebase '/h_iter'], 'epf')
    end
end