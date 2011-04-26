pathfile = './';
dark = readtimeseries([pathfile 'R_20110413/Dark_1/img_000000*__000.tif']);
bright = readtimeseries([pathfile 'R_20110413/Bright_1/img_000000*__000.tif']);
sv= size(bright);
out_cal = cal_readnoise(bright, dark(0:sv(1)-1,0:sv(2)-1,:));