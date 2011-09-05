sep = .4;
f1_2D=makeGauss2D(x,0,p.sig1);
f2_2D=makeGauss2D(x,sep,p.sig2);
int1=500;
bgvec = [0,100,500,1000];
% i500=zeros(10,11,3);
% for ii=1:length(bgvec)
%     o=bgvec(ii);i500(:,:,ii)=double(noise(resample(int1*f1_2D+int1*f2_2D+o,.2), 'poisson'));
% end
% 
% if savethis
%     for ii=1:3
%         dipshow(i500(:,:,ii));
%         name = ['twoSources_Sep1_int' num2str(int1) '_bg' num2str(bgvec(ii)) '_sep' num2str(sep)];
%         SaveImageFULL(['images/' name]);
%     end
% end

int1=1000;
i1000=zeros(10,11,3);
for ii=1:length(bgvec)
    o=bgvec(ii);i1000(:,:,ii)=double(noise(resample(int1*f1_2D+int1*f2_2D+o,.2), 'poisson'));
end


if savethis
    for ii=1:4
        dipshow(i1000(:,:,ii));
        name = ['twoSources_Sep1_int' num2str(int1) '_bg' num2str(bgvec(ii)) '_sep' num2str(10*sep)];
        SaveImageFULL(['images/' name]);
    end
end