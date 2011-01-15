function yfit=gaussian_merge(beta,xdummy)
  global xpix ypix wbox;
 
z=beta(3)*exp(-2*((xpix-beta(1)).*(xpix-beta(1))+(ypix-beta(2)).*(ypix-beta(2)))/(beta(4))^2);
for i=1:wbox
    for j=1:wbox
      k=(i-1)*wbox+j;
      yfit(k)=z(i,j);
    end
end

