function I=fisherInformationMatrix(xhires,f1,f2, int1, int2, pixelversion, offset,osf)
% =fisherInformationMatrix(xhires,f1,f2, int1, int2, pixelversion,offset,osf)

if ~exist('pixelversion','var')
    pixelversion = 0; 
    osf=1;
end
if ~exist('offset','var')
    offset=0;
end

if ndims(xhires)==2 %1D vector
    dx=xhires(2)-xhires(1);
else
    dx=xhires(1,2,1)-xhires(1,1,1);
end

if pixelversion
    % PSFs:
    f1_2d=reshape(f1,size(xhires,1),size(xhires,2)); 
    f2_2d=reshape(f2,size(xhires,1),size(xhires,2)); 
    fp1_2d= binsumImage(f1_2d,[osf,osf]);
    fp2_2d= binsumImage(f2_2d,[osf,osf]);
    fp1=fp1_2d(:);
    fp2=fp2_2d(:);
    % derivatives:
    gfp1_2d= binsumImage(gradient(f1_2d,dx),[osf,osf]);
    gfp2_2d= binsumImage(gradient(f2_2d,dx),[osf,osf]);
    gfp1=gfp1_2d(:); 
    gfp2=gfp2_2d(:); 
else     
    fp1=f1; 
    fp2=f2;
    gfp1=gradient(f1,dx);
    gfp2=gradient(f2,dx);
%     xp=x;
end
I=zeros(2); %Fisher Information matrix 2x2;
p1 = int1*fp1+int2*fp2+offset;
% p2 = int2*fp1+int1*fp2+offset;

prec = eps;
mask = abs(gfp1.^2)>prec;
kernel11 = (1./p1(mask)).*(int1*gfp1(mask)).^2;
mask = abs(gfp2.^2)>prec;
% kernel22 = (1./p2(mask)).*(int2*gfp2(mask)).^2;
kernel22 = (1./p1(mask)).*(int2*gfp2(mask)).^2;
mask = abs(gfp1.*gfp2)>prec;
kernel12 = (1./p1(mask))*int1.*int2.*gfp1(mask).*gfp2(mask);
  
I(1,1)=sum(kernel11);
I(2,2)=sum(kernel22);
I(1,2)=sum(kernel12);
I(2,1)=I(1,2);