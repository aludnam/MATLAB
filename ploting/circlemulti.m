function circlemulti(centers,radius,nop,style)

% H=CIRCLEMULTI(CENTER,RADIUS,NOP,STYLE)
% Draws multi circles specified by centres as rows of the matrix centres.
% example: circlemulti(10*rand(10,2), 3,20,'--r')
hold on
sc=size(centers,1);
if length(radius)==1
    radius = repmat(radius,sc,1);
end
for ii=1:sc
    circle(centers(ii,:),radius(ii),nop,style);
end