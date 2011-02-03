function [dpixc_n, dveccr_n, blinkmat_n, p_n] = catsimul(namedir, catitervec)
dpixc_n=[];
dveccr_n=[];
blinkmat_n=[];
p_n=[];

for ii=catitervec
    namefile = [namedir '-iter_' num2str(ii)];
    load ([namefile '.mat'])
    dpixc_n=cat(3,dpixc_n, dpixc);
    dveccr_n=cat(2,dveccr_n, dveccr);
    blinkmat_n=cat(2,blinkmat_n, blinkmat);
    
end
p_n=p;
p_n.Nt=p.Nt*length(catitervec);




