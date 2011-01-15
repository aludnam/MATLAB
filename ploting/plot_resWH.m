showbg = 0;
showpoints = 1;
ca
mh = mean (res.h,2);
[mhs,ihs] = sort(mh,'descend');

wr = reshape(res.w,peval.nx,peval.ny,size(res.w,2));
mkdir ('w')

sw = size(res.w,2);
a = ceil(sqrt(sw));
b = ceil(sw/a);
for ii=1:sw
    dipshow(ii,wr(:,:,ihs(ii)))
    SaveImageFULL(['w/w_' num2str(ii)])
    figure(sw+1) 
    subplot(a,b,ii)
    imagesc(wr(:,:,ihs(ii)))
    set(gca,'dataaspectratio',[1 1 1])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    xlabel(num2str(ii))
end

if showpoints
    for ii=1:size(res.w,2)-1
        dipshow(ii,wr(:,:,ihs(ii+1)))
        hold on
        scatter(p.x_vec-.5, p.y_vec-.5,'r')
        hold off
        SaveImageFULL(['w/w2_' num2str(ii)])
    end
end

figure
imagesc(res.h(ihs(2:end,:),:))
if ~showbg; ylim([0.5 peval.ncomp+0.5]); end
xlabel('slice # (time)')
xlabel('component')
SaveImageFULL(['Himage'])

hsort=res.h(ihs,:);

cc = corrcoef(hsort(2:end,:)');
figure
imagesc(cc)
set(gca, 'DataAspectRatio',[1 1 1])
colorbar
SaveImageFULL('Hcorrel')

figure
imagesc(abs(cc))
set(gca, 'DataAspectRatio',[1 1 1])
colorbar
SaveImageFULL('Hcorrelabs')

figure
bar(mhs(2:end))
xlabel('component - last is background')
ylabel('mean (H,2)')    
grid on
% if ~showbg; xlim([1.5 peval.ncomp+1.5]); end
SaveImageFULL(['Hmean_bar'])
