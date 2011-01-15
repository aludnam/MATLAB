% aaa= readtimeseries( '/Volumes/O_MANDULA/data/qdots/R1/100222/7W5UG85Y_F00000011.tif');
aaa= readtimeseries( '/Users/ondrejmandula/Documents/data/100222/7W5UG85Y_F00000011.tif');
aaa14=squeeze(aaa(:,:,0:100,0));
aaa14small=aaa14(x_offset:x_offset+x_size, y_offset:y_offset+y_size,:);
aaad=double(aaa14small);
imagesc(sum(aaad,3))
colormap('gray')
hold on
scatter(xf_all-1, yf_all-1,'xr')
hold off