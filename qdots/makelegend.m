function makelegend (offset, namedisplay, xlabel_my, ylabel_my)
% function makelegend (offset, namedisplay, xlabel_my, ylabel_my)
% sep = .1:.1:2;
% offset = [10 100 1000];
kk = 1;
% t{1} = namedisplay;

% t{2} = 'ICA';

for ii=1:1
    for jj= 1:length(offset)
%         l{kk} = [t{ii} 'offset ' num2str(offset(jj))];
%         l{kk} = [t{ii} num2str(offset(jj,:))];
%         l{kk}=namedisplay{jj};
        l{kk}=[namedisplay  num2str(offset(jj))];
        kk = kk+1;
    end
end

legend(l)
xlabel (xlabel_my)
ylabel (ylabel_my)
% xlabel('separation [pixels]')
% ylabel('loc err [pixels]')
grid on
% xlim([0 10])
% ylabel('loc err / separation')