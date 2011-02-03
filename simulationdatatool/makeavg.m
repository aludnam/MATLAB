function dpixc = makeavg(dpixc_ind, blinkmat, offset, avg)
% dpixc = makeavg(dpixc_ind, blinkmat, offset, avg)      
% avg = 0 -> noise free image
if size(blinkmat,2) == 1 %one slice only
    dpixc_nonoise = dpixc_ind'*blinkmat + offset;
else
    dpixc_nonoise = array2im(dpixc_ind'*blinkmat + offset);
end
dpixc_dip = newim(dpixc_nonoise);
for ii=1:avg
    dpixc_dip = dpixc_dip + noise(dpixc_nonoise,'poisson');    
end
if avg==0
    dpixc_dip = dpixc_nonoise;
end
navg=max(1,avg);
dpixc = double(dpixc_dip/navg);

