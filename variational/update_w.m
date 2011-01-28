function w=update_w(dvevc, n)
% w=update_w(dvevc, n)
% Computes update for w_jk in variational approximation (Buntine & Jakulin
% DCA 2006)
% dvec: data - each column is one image (#pixels=J x #documents=I)
% n: (J x K x I )

i=size(dvec,2); % number of images


