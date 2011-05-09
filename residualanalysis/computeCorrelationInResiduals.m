function out = computeCorrelationInResiduals(data_orig, data_eval)


resid = data_orig - data_eval;
ccd = corrcoef(resid');
out = max(max(ccd - eye(size(ccd))));
% 
% 
% method='average';
% data=resid_norm;
%         sized = size(data);
%         dveccr= reshape(data,sized(1)*sized(2), sized(3));
%         ccd = (corrcoef(dveccr'));
%         ccds = correlation2distance(ccd, 'pos');
%         
%         Z = linkage(ccds,method);
%         
%         if savethis
%             save ([res_dir '/corrcoefdist.mat'], 'Z')
%         end
%         maxcorrel(ii,iteration)=1-min(Z(:,3));