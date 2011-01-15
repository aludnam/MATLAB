sig = 1.20644;
offset = 64.4208;
maxphot = 1.4027e+03;
g=maxphot*dip_image(makegauss([50 50], sig, [100 100]));
gn=noise(g+offset,'poisson');


offset = 92.8808;
sig = 1.30246;
maxphot = 2580.42;

g=maxphot*dip_image(makegauss([50 50], sig, [100 100]));
gn=noise(g+offset,'poisson');

