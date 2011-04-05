for k=[2 3 28 29]
    
    [meanf,varf,entf, kurt]=computeMeanVarFromPermutations(f',k);
    %         fpix=squeeze(oroi3(indexpixel,:,:));
    %         [meanf,varf,entf]=computeMeanVarFromPermutations(fpix,k);
    [color, marker, line] = getcolorfromindex(gca, k);
    %         figure(indexpixel)
    hold on
    %         scatter(meanf,sqrt(varf),20,color)
    scatter(meanf,kurt,20,color)
    l{k}=['k=' num2str(k)];
    hold off
    
end