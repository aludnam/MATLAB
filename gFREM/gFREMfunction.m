function y=gFREMfunction(x,f1,f2, int1, int2, showim)
% y=gFREMfunction(x,f1,f2, int1, int2, showim)

if ~exist('showim','var')
    showim = 0;
end

p = int1*f1+int2*f2;
dx=x(2)-x(1);
kernel = (1./p).*(int1*gradient(f1,dx)-int2*gradient(f2,dx)).^2;
y=trapz(x,kernel);

if showim
figure; 
plot(x,int1*f1, x,int2*f2,x,p,':');
hold on
plot(x,kernel,'r--');
grid on
end




