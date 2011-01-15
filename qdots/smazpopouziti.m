[dpixc, dpixc_ind, blinkmat, peval, p] = readSimulMulti(peval.N, 100, 1, 1, peval);
peval.dir_res = [peval.path_res peval.namedir_data];
bicp = computebic(peval.dir_res, [peval.namedir_data '_iter1_avg1_nc'], ncomp_vec, dpixc);
plot(ncomp_vec, bicp);
xlabel('number of components');
ylabel('BIC score');
grid on; 
xlim([5,30])
SaveImageFULL('BIC')
