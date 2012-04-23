a=readim('img_000000000__000.tif')

ar=a(293:415,35:155)

allshiftsX=zeros(1,1118)
allshiftsY=zeros(1,1118)

maxN=1117;
myStep=10;

for n=1:myStep:maxN
    b=readim(sprintf('img_00000%.4d__000.tif',n));
    br=b(293:415,35:155);
    res=kcorrelator({ar,br},{'v','Myshift.txt'});
    myshift=load('myShift.txt');
    allshiftsX(n)=myshift(1);
    allshiftsY(n)=myshift(2);
    fprintf('Processing: %d\n',n);
end
pixelsize=16000/157.5;
plot(pixelsize*allshiftsX(1:myStep:maxN))
hold on
plot(pixelsize*allshiftsY(1:myStep:maxN))
xlabel('Frame #')
ylabel('Drift [nm]')
