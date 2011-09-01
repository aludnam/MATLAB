function I=fisherInformationMatrix(x,f1,f2, int1, int2, pixelversion, offset)
% I=fisherInformationMatrix(x,f1,f2, int1, int2, pixelversion, offset)

if ~exist('pixelversion','var')
    pixelversion = 0; 
end
if ~exist('offset','var')
    offset=0;
end

if ndims(x)==2 %1D vector
    dx=x(2)-x(1);
else
    dx=x(1,2,1)-x(1,1,1);
end

if pixelversion    
    [fp1, xp]= pixelize1d(x,f1);
    fp2= pixelize1d(x,f2);
    gfp1=pixelize1d(x,gradient(f1,dx));
    gfp2=pixelize1d(x,gradient(f2,dx));
else     
    fp1=f1; 
    fp2=f2;
    gfp1=gradient(f1,dx);
    gfp2=gradient(f2,dx);
    xp=x;
end
I=zeros(2); %Fisher Information matrix 2x2;
p1 = int1*fp1+int2*fp2+offset;
p2 = int2*fp1+int1*fp2+offset;

prec = 10^-8;
mask = abs(gfp1.^2)>prec;
kernel1 = (1./p1(mask)).*(int1*gfp1(mask)).^2;
mask = abs(gfp2.^2)>prec;
kernel2 = (1./p2(mask)).*(int2*gfp2(mask)).^2;
mask = abs(gfp1.*gfp2)>prec;
kernel3 = (1./p1(mask))*int1.*int2.*gfp1(mask).*gfp2(mask);
  
%%% test
I(1,1)=sum(kernel1);
    I(2,2)=sum(kernel2);
    I(1,2)=sum(kernel3);
    I(2,1)=I(1,2);
%%% 

% if pixelversion
%     I(1,1)=sum(kernel1);
%     I(2,2)=sum(kernel2);
%     I(1,2)=sum(kernel3);
%     I(2,1)=I(1,2);
% else
%     
%     I(1,1)=trapz(xp(1:length(kernel1)),kernel1);
%     I(2,2)=trapz(xp(1:length(kernel2)),kernel2);
%     I(1,2)=trapz(xp(1:length(kernel3)),kernel3);
%     I(2,1)=I(1,2);
% end