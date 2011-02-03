function blinkmat = blinkmat_markov(N,Nt, intensity_vec, probtrans)
intensity_mat = repmat(intensity_vec, 1,Nt);
changemat = rand(N,Nt)<probtrans;
statemat = mod(cumsum(changemat,2),2);
initvec = rand(N,1)>0.5; %initial state of the blinkmat
ivt = ~(initvec == statemat(:,1));
statematinit = mod(statemat+repmat(ivt,1,Nt),2);
blinkmat = statematinit .* intensity_mat; %different intensities...
