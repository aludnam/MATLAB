
% int1_cell2mat=reshape(cell2mat(int1_mat),length(int1_vec),length(int1_multi));
% pint1_cell2mat = reshape(cell2mat(pint1_mat),length(int1_vec),length(int1_multi));

figure; 
hold on
co = get(gca, 'colororder');

mm= max(cellfun('length', int1_mat));
im1=zeros(mm,length(int1_mat));
pi1=zeros(mm,length(int1_mat));

for ii=1:length(int1_mat)
    im1(1:length(int1_mat{ii}),ii)=int1_mat{ii};
    pi1(1:length(pint1_mat{ii}),ii)=pint1_mat{ii};
%    plot(int1_mat{ii},pint1_mat{ii},'-o','color',co(ii,:))
%     bh=bar(int1_mat{ii},pint1_mat{ii},.1),%,.1/length(int1_mat{ii}));
%     set(bh,'facecolor',co(ii,:))
end
h=plot(im1,pi1,'o-');
grid on
xlim([min(im1(:))-.2,max(im1(:))+.2])
hold off
xlabel('intensity')
ylabel('probability')
grid on
legend(l)
setfontsizefigure(12)
setforsave(gcf,1.1)