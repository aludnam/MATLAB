function blinkmat = blinkmatpowerlaw(ncomp, nt, alpha, intensities)
% blinkmat = blinkmatpowerlaw(ncomp, nt, alpha, intensities)
% Generates the binkmat distributed according to the power law
% blinkmat: ncomp x nt matrix of intensities
% ncom: number of components
% nt: number of time slices
% alpha: power law coefficient (close to 1.5-1.7)
% intensities: vector of the individual intensities of the individal
% sources

if ~exist('intensities', 'var')
    intensities = ones(1,ncomp);
end

% alpha = 1.7;
Ntimes = 10^6;
rx=rand(1,Ntimes);
tx=ceil((alpha-1)*rx.^(1-alpha));
% ncomp = 10;
nt=1000;

blinkmat=zeros(nt,ncomp); %transposed for computation convenience
value=[0 1];
indexstart=1;
indexstop=1;
jj=1;

while indexstart+tx(jj) <= numel(blinkmat)

indexstop=indexstart+tx(jj)-1;
blinkmat(indexstart:indexstop)=value(mod(jj,2)+1);
indexstart=indexstop+1;
jj=jj+1;
end


blinkmat=bsxfun(@times, makeveccolumn(intensities), blinkmat');

imagesc(blinkmat)