sep =[.1 .2 .5];
offset = [100 1000];
niter = 1;
avg_num = [0 1 2 10 50]; %number of averages os samples (reducing Poisson noise)
% avg_num=0 --> noise free image...

path_res = [cd '/'];
path_data = '~/project/data/qdots/S44/';
prename = 'S44_sep_';

cd (path_data)
savethis = 1;
separ_Sexample(sep, offset, path_data, path_res, prename, niter, avg_num, savethis);