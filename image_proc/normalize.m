function out = normalize(in,p)
% out = normalize(in,p)
% normalizes 1d vector or 2d image or each slice of 3D image in the p-norm
% (sum(in.^p) = 1; (p=1 by default)
% for p==0 the image is normalized to the maximum value==1
isdip=0;
if strcmpi(class(in),'dip_image')
    isdip=1;
    in=double(in); 
end

if nargin < 2 
    p=1;
end
nd=ndims(in);
if p==0
    fprintf('Normalizing to maximum.\n')
else 
    fprintf('Normalizing to L%d norm.\n',p)
end
    
switch nd
    case 1 %1D        
        if p==0
            out = in/max(in);           
        else
            out = in/sum(in.^p);                
        end
    case 2 %2D
        if p==0
            out=bsxfun(@rdivide, in,max(in(:),[],1));
        else
            out=bsxfun(@rdivide, in,sum(in(:).^p,1));
        end
    case 3 %3D        
        si = size(in);
        inr = reshape(in,si(1)*si(2),si(3));
        if p==0
            outr = bsxfun(@rdivide, inr, max(inr,[],1));            
        else
            outr = bsxfun(@rdivide, inr, sum(inr.^p,1));            
        end
        out = reshape(outr, si);
end
if isdip
    out=dip_image(out); 
end
