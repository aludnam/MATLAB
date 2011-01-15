% out2 = readtimeseries('/afs/inf.ed.ac.uk/user/s08/s0880377/project/data/qdots/R1/100222/7W5UG85Y_F00000014.tif','',[],1,1);
% out14=squeeze(out2(:,:,:,0));
% o=dip_image(out14);
o=dip_image(data_out);
box1=4;
box2=2;
nclick = 10;
dipshow(7,mean(o,[],3));
coords = dipgetcoords(7,nclick);
box1x1 = coords(:,1) - box1;
box1x2 = coords(:,1) + box1;
box1y1 = coords(:,2) - box1;
box1y2 = coords(:,2) + box1;
for ii=1: 7%length(coords)
    oroi1 = o(box1x1(ii):box1x2(ii), box1y1(ii):box1y2(ii),:);
    [maxval, maxcoord] = max(max(oroi1,[],3));
    box2x1 = maxcoord(1)-box2;
    box2x2 = maxcoord(1)+box2;
    box2y1 = maxcoord(2)-box2;
    box2y2 = maxcoord(2)+box2;
    oroi2= oroi1(box2x1:box2x2,box2y1:box2y2,:);
    f(ii,:) = squeeze(double(sum(oroi2,[],[1 2])));
    
end

index=1;
figure(1)
for ii=1:N    
    subplot(N,2,index)
    plot(f(ii,:))
      if ii<N 
        set(gca,'xtick',[],'ytick',[])
    else 
        set(gca,'ytick',[])
    end
    index = index +1;
    subplot(N,2,index)
    hist(f(ii,:),50);
    xlim([1.5,4]*10^5)
    if ii<N 
        set(gca,'xtick',[],'ytick',[])
    else 
        set(gca,'ytick',[])
    end
    index = index +1;
end
SaveImageFULL(['block' num2str(q) '_timeprofilesingledot'])
figure(2)
bar(max(f,[],2),'r')
hold on
bar(mean(f,2),'g')
bar(min(f,[],2),'b')
errorbar([1:8], mean(f,2), sqrt(var(f,[],2)),'+k')
l{1}='max';
l{2}='mean';
l{3}='min';
l{4}='sqrt(var)';
legend(l);
grid on
hold off
SaveImageFULL(['block' num2str(q) '_statistics'])