function [vard, I3d]=computeSeparationVariance(x,l1,l2,sig,int_vec, pint, pixelizeversion, offset)
% [vard, I3d]=computeSeparationVariance(x,l1,l2,sig,int_vec, pint, pixelizeversion)
% Computes variance (as a inverse of the Fisher Information) of two points
% separated by a distance d=l1-l2.

if ~exist('offset','var')
    offset = 0; % background 
end

int1_vec=int_vec(1,:);
int2_vec=int_vec(2,:);
pint1=pint(1,:);
pint2=pint(2,:);


lint1=length(int1_vec);
lint2=length(int2_vec);
ll = length(l2);
vard = zeros(1, ll);
I3d = zeros(2,2, ll);

f1=makeGauss(x,l1,sig(1));                  % creates PSF (gauss approx)
for ind_dist=1:ll                           % distance
    I=zeros(2);    
    f2=makeGauss(x,l2(ind_dist),sig(2));    % creates PSF (gauss approx) shifted to l2
    for ind_int1=1:lint1                    % intensity of the source 1
        int1=int1_vec(ind_int1);
        for ind_int2=1:lint2                % intensity of the source 1
            int2=int2_vec(ind_int2);            
            % This is accumulating the Fisher Information matrix for
            % different intensities:
            if ~and(int1==0,int2==0)
                vard(ind_dist)=test_fisherInformationMatrix(x,f1,f2, int1, int2,pixelizeversion, offset);            
            end
        end
    end
%     I3d(:,:,ind_dist)=I;    
%     vard(ind_dist)=[1,-1]/I*[1,-1]';    
I3d=0;

end
q=0;
