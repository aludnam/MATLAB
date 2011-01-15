function varargout=testNMFS51(varargin)
% function [X1,X2,dh,minX]=testNMF2
% function [X1,X2,dh,minX]=testNMF2(w,h,p)
v=varargin{1};
w=varargin{2};
h=varargin{3};
hinit=varargin{4};
h_dovec=varargin{5};
p=varargin{end};

[n,m]=size(w*h);
npix = p.nx*p.ny;
wbg = (1/npix)*ones(npix,1);
p.maxh=100;
checkvec = [10, 100, 500 10000 p.maxh]; %where the values are recorded

if max(p.separ) < 0.5
    [X1, X2] = meshgrid(14.8:0.05:15.5, 14.8:0.05:15.5);
else
    [X1, X2] = meshgrid(14.8:0.1:16, 14.8:0.1:16);
end

for ixmat=1:size(X1,1)
    for jxmat=1:size(X1,2)
        x=[X1(ixmat, jxmat), 15, X2(ixmat, jxmat), 15];
        wg = makegauss(x, p.s, [p.nx p.ny]);
        w = reshape(wg, p.nx*p.ny,size(wg,3));
        % normalization of all w:
        sumw = sum(w,1);
        w = w./repmat(sumw,n,1); %normalization of each component
        h(1:2,:) = h(1:2,:).*repmat(sumw',1,m); %to keep the multiplication equal
        
        
        options = zeros(1,18);
        options (1)=1; %to display error values
        options (7)=1;
        options(9)=0; %to check gradient
        options(14)=p.maxh; %maximum number of iterations
        he=log(hinit(1:2,:));
        [rHkte, options, flog, pointlog] = conjgrad('ddivHexpS51', reshape(he,1,m*2), options, 'gradddivHexpS51', v, w, wbg, hinit(3,:));
        %S52:[rHkte, options, flog, pointlog] = graddesc('ddivHexpS51', reshape(he,1,m*2), options, 'gradddivHexpS51', v, w, wbg, hinit(3,:));
        dh(ixmat,jxmat)=ddivergence(v, [w,wbg]*[exp(rHkte');hinit(3,:)]);
    fprintf('[%g\t%g]\tDdivergence %g\n',X1(ixmat,jxmat), X2(ixmat,jxmat), flog(end))
    end
    %         mhd(ixmat,jxmat,:,:)=mean(htr,3); %4dimension: [x,y,z,components]
    
end
    
% for ii=1:zxmat-1
%     imtmp = dh(:,:,ii);
%     miXval(ii) = min(imtmp(:));
%     [mx(ii),my(ii)]=find(miXval(ii)==imtmp,1,'first');
%     minX(ii,:)=[X1(mx(ii),my(ii)), X2(mx(ii),my(ii))]
% end
varargout{1}=X1;
varargout{2}=X2;
varargout{3}=flog;
varargout{4}=[];
varargout{5}=[];
varargout{6}=dh;
varargout{7}=[];
varargout{8}=[];
varargout{9}=p;


% function varargout=assignvalues(varargin)
% v=varargin{1};
% w=varargin{2};
% h=varargin{3};
% ixmat=varargin{4};
% jxmat=varargin{5};
% zxmat=varargin{6};
%
% dh(ixmat,jxmat,zxmat)=ddivergence(v, w*h);
% htr(zxmat,:,:) = h;
% varargout{1}=dh;
% varargout{2}=htr;
% end