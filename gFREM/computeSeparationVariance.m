function vard=computeSeparationVariance(x,l1,l2,sig,int_vec, pint, pixelizeversion)
% vard=computeSeparationVariance(x,l1,l2,sig,int_vec, pint, pixelizeversion)
% Computes variance (as a inverse of the Fisher Information) of two points
% separated by a distance d=c1-c2.

sig1=sig(1); sig2=sig(2);
int1_vec=int_vec(1,:);
int2_vec=int_vec(2,:);
pint1=pint(1,:);
pint2=pint(2,:);

lint1=length(int1_vec);
lint2=length(int2_vec);

vard = zeros(1, length(l2));
for ind_dist=1:length(l2)           % distance
    I=zeros(2);
    for ind_int1=1:lint1            % intensity of the source 1
        int1=int1_vec(ind_int1);
        for ind_int2=1:lint2        % intensity of the source 1
            int2=int2_vec(ind_int2);
            f1=makeGauss(x,l1,sig1);            % creates PSF (gauss approx)
            f2=makeGauss(x,l2(ind_dist),sig2);  % creates PSF (gauss approx) shifted to l2
            % This is accumulating the Fisher Information matrix for
            % different intensities.
            if ~and(int1==0,int2==0)
                I=I+fisherInformationMatrix(x,f1,f2, int1, int2,pixelizeversion);            
            end
        end
    end
    I3d(:,:,ind_dist)=I;
    vard(ind_dist)=[1,-1]/I*[1,-1]';
end
q=0;
