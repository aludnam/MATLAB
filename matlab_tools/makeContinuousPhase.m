% Algorithm to extract the phase information and add a global phase to avoid phase wraps (unwrapping if possible).
% function myangle=makeContinuousPhase(PF,mask);
% PF = amplitude distribution
% myangle = unwrapped phase map
% Author: Rainer Heintzmann
function myangle=makeContinuousPhase(PF,mask);
if nargin < 2
    myangle=angle(PF);   
    PF = PF;
else
    myangle=angle(PF).*mask;    %put mask to ignore all around the aperture
    PF = PF .* mask;
end
mydx = dx(real(PF)) + i* dx(imag(PF));  % derive complex with respect to X
mydy = dy(real(PF)) + i* dy(imag(PF));  % derive complex with respect to Y

x_grad_phi = mydx ./ (i*PF);   % convert to gradient over phi instead of PF
y_grad_phi = mydy ./ (i*PF);

x_grad_phi(~mask)=0;
y_grad_phi(~mask)=0;

mykx = xx(size(PF),'freq')*2*pi;    %create image of same size as PF, 
myky = yy(size(PF),'freq')*2*pi;    %unwraped image?

mask_x = (mykx ~= 0);       %make mask to cover infinity gap of kx and ky
mask_y = (myky ~= 0);

divisor = mykx*mykx+myky*myky;      %for simplification

FT_x_grad_phi=ft(x_grad_phi);       %make fourier image
FT_y_grad_phi=ft(y_grad_phi);

integrated = dip_image(PF*0,'scomplex');   
% intrgration over x direction
integrated(mask_x) = mykx(mask_x) .* mykx(mask_x) .* FT_x_grad_phi(mask_x) ./ (i.*mykx(mask_x).*divisor(mask_x));
% integration over y direction
integrated(mask_y) = integrated(mask_y) + myky(mask_y) .* myky(mask_y) .* FT_y_grad_phi(mask_y) ./ (i.*myky(mask_y).*divisor(mask_y));
% take only real numbers, imaginary are only small mistakes
integrated = real(ift(integrated));


% Take the middle point as a reference to determine the additive constant
integrated = integrated + myangle(floor(size(PF,1)/2),floor(size(PF,2)/2)) - integrated(floor(size(PF,1)/2),floor(size(PF,2)/2));
%
mydiff=myangle - integrated;
mydiffclear=2*pi*floor((mydiff+pi)/2/pi);

myangle=myangle-mydiffclear; % and this is the result

    