function V = cbcldata
% cbcldata - read face image data from cbcl database
%

global imloadfunc;
    
% This is where the cbcl face images reside
thepath = '../data/cbcl-face-database/face/';

% Create the data matrix
V = zeros(19*19,2429);

% Read the directory listing
D = dir(thepath);

% Step through each image, reading it into the data matrix
% Note: The (+2) is just to avoid '.' and '..' entries
fprintf('Reading in the images...\n');
for i=1:2429,
    switch imloadfunc,
     case 'pgma_read',
      I = pgma_read([thepath D(i+2).name]);
     otherwise,
      I = imread([thepath D(i+2).name]);
    end
    V(:,i) = reshape(I,[19*19, 1]);
    if rem(i,100)==1, fprintf('[%d/24]',floor(i/100)); end
end

fprintf('\n');

% Same preprocessing as Lee and Seung
V = V - mean(V(:));
V = V / sqrt(mean(V(:).^2));
V = V + 0.25;
V = V * 0.25;
V = min(V,1);
V = max(V,0);


% Additionally, this is required to avoid having any exact zeros:
% (divergence objective cannot handle them!)
V = max(V,1e-4);

