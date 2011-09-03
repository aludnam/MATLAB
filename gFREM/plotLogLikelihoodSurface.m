bg = 100;
% dvec=[0 1 2 3 4 5];
dvec = 1;
fontsize = 15; 

for ii=1:length(dvec)
    p.d = dvec(ii);
    
    [logl, c1, c2, p] = computeLogLikelihood(p.d,bg);
    
    figure;
    % [DX,DY]=gradient(logl, c1(2)-c1(1));
    % quiver(c1,c2,DX,DY)
    % hold on
    contour(c1,c2,logl, 20)
    hold on
    scatter(p.c1_true, p.c2_true,'*k');
    scatter(p.c2_true, p.c1_true,'*k')
    xlabel('c_1 (center of the source one)')
    ylabel('c_2 (center of the source two)')
    colorbar
    set(gca,'dataaspectratio',[1 1 1])
    grid on
    setfontsizefigure(fontsize)
    % colormap(gray)
    if savethis
        name = ['images/LogLikelihoodSurface_d' num2str(p.d) '_bg' num2str(p.bg)];
        SaveImageFULL(name)
    end
    
    
    figure;
    mesh(c1,c2,logl)
    hold on
    scatter3(p.c1_true, p.c2_true, max(logl(:)),'*k');
    scatter3(p.c2_true, p.c1_true, max(logl(:)),'*k')
    % title(['bg=' num2str(p.bg)]);
    
    xlabel('c_2')
    ylabel('c_1')
    zlabel('log likelihood')
%     zlim([-2.5e4 0])
    setfontsizefigure(fontsize)
    % colormap(gray)
    if savethis > 2
        name = ['images/LogLikelihoodSurface3d_d' num2str(p.d) '_bg' num2str(p.bg)];
        SaveImageFULL(name)
    end
    
    % figure;
    llcut(:,ii)=(logl(find(fliplr(eye(size(logl))))))';
    % mm=mm+1;
end