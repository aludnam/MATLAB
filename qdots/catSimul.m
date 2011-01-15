function [dpix_n, blinkmat_n] = catSimul(path_data, namedir_data, catitervec)
dpix_n=[];
blinkmat_n=[];

for ii=catitervec
    namefile = [namedir_data '-iter_' num2str(ii)];
    load ([path_data namedir_data '/' namefile])
    dpix_n=cat(3,dpix_n, dpixc);
    blinkmat_n=cat(2,blinkmat_n, blinkmat);    
end



