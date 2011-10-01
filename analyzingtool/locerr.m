function [dist, indexmutual]=locerr(xt,yt,xl,yl,showimage)
% [dist, indexmutual]=locerr(xt,yt,xl,yl,showimage)
% Estimates minimum distances between the localized points (xl,yl) and the
% true locations (xt,yt) in the vector dist
if ~exist('showimage','var')
    showimage = 0; 
end
sl=size(xl);
if sl(1)<sl(2) %making column vectors
    xl=xl'; yl=yl';
end
st=size(xt);
if st(1)<st(2) %making column vectors
    xt=xt'; yt=yt';
end
st=length(xt);
mt=[xt,yt];
ml=[xl,yl];
m=[mt;ml];
mdist=squareform(pdist(m)); %distance matrix
mdistmutual=mdist(st+1:end, 1:st); %mutual distances
[dist, indexmutual]=findminimaldist(mdistmutual);
if showimage
    scatterpoints(xt,yt,xl,yl,indexmutual)
end
end

function [md, rowsort] = findminimaldist(A)
%find minimum mutual distances and makes sure that each point is allocated
%only one localized position
Aorig=A;
sv = size(A);
for ii=1:min(sv)
    [md(ii), index] = min(A(:));
    [minrow, mincol] = ind2sub(size(A),index);
    A = removerows(A,minrow);         
    A = removerows(A',mincol);
    A=A';
    [irow_tmp, icol_tmp]=find(Aorig==md(ii));
    irow(ii)=irow_tmp(1);
    icol(ii)=icol_tmp(1);
end
[Y,I]=sort(icol);
rowsort=irow(I);
end

function scatterpoints(xt,yt,xl,yl,indexmutual)
lindex = length(indexmutual);
figure; hold on
co=get(gca, 'colororder');
lco=length(co);
rat = lindex/lco;

if rat > 1 %more clusters than colors
    co = repmat(co, ceil(rat), 1);
end

for ii=1:lindex
    scatter(xt(ii), yt(ii),[],co(ii,:));
    scatter(xl(indexmutual(ii)),yl(indexmutual(ii)),[],co(ii,:),'x');
    line([xt(ii), xl(indexmutual(ii))], [yt(ii),yl(indexmutual(ii))],'color',co(ii,:));
    grid on
    l{1}='true'; l{2}='estimate';
    legend(l)
end
end
