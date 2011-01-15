function plotHexptiled(res,holdonyes,col)
if nargin<2
    holdonyes=0;
end
if nargin<3
    col=[0 0 1];
end
htrue_val = res{1}.p.meanblinkmat*res{1}.p.htrue;
margin = [1000 1000];
for ll=1:4
    subplot(2,2,ll)
    if holdonyes
        hold on
    end
    nsteps=length(res{ll}.htrace);
    plot(repmat([0:nsteps-1]',1,2), exp(res{ll}.htrace),'Color',col)
    hold on
    plot(repmat([1,res{ll}.p.maxh-1]',1,2), repmat(res{ll}.p.meanblinkmat*res{ll}.p.htrue,1,2)','--r');
    hold off
    
    xlabel('iter');
    ylabel('h')
    ylim([min(htrue_val)-margin(1), max(htrue_val)+margin(2)])
    set(gca, 'XScale','log','XMinorTick','on','XMinorGrid','on')
%     set(gca, 'XTick',[0,10,100,1000, 10000, 100000])
    grid on
end
hold off