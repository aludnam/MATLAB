function plot_log_normal( x,params,hAx,plot_num,fontsize )
% plot the log-normal distribution with parameters "m" and "s"
% 
% the distribution is given by:
%
%        p(x) = sqrt(1/(2*pi))/(s*x)*exp(- (log(x-m)^2)/(2*s^2))
%
% format:   plot_log_normal( x,params,hAx,plot_num,fontsize )
%
% input:    x         - X axis, for the plot
%           params    - the distribution parameter, RMS error and VAR or CRB
%           hAx       - where to plot the distribution curve
%           plot_num  - since the curve is added with a text to the axes,
%                       this parameter specifies where the text should be displayed
%                       and what color to choose for the curve
%           fontsize  - size of the font of the text, default 9
%
%
% example:  plot_log_normal( x,fit_ML_log_normal(data),hAx,3 )
%

% init graphic parameters
switch plot_num
case 1, cl = [1 0 0];
case 2, cl = [0 1 0];
case 3, cl = [1 0 1];
case 4, cl = [0 1 1];
case 5, cl = [0.5 0.5 0];
end
if ~exist('fontsize')
    fontsize = 9;
end
p       = plot_num*0.15 + 0.3;
ylimit  = ylim(hAx);
xlimit  = xlim(hAx);
fnc_txt = '1/(\bfs\cdotx\rm\surd\it{2\pi}) \bf{\cdot e}^{-(ln \bf{x - m}\rm)^2/2\bf{s^2}}\rm';
if isfield( params,'VAR' )
    txt     = sprintf( '\\fontsize{%d}\\bfLog-Normal PDF\\rm with %s:  %s\nm = %g    s^2 = %g\nVAR(mu) = %1.3g  VAR(s^2) = %1.3g\nRMS err = %1.3g\n',...
        fontsize,params.type,fnc_txt,params.m,params.s2,params.VAR_m,params.VAR_s2,params.RMS );
else
    txt     = sprintf( '\\fontsize{%d}\\bfLog-Normal PDF\\rm with %s:  %s\nm = %g    s^2 = %g\nCRB(mu) = %1.3g  CRB(s^2) = %1.3g\nRMS err = %1.3g\n',...
        fontsize,params.type,fnc_txt,params.m,params.s2,params.CRB_m,params.CRB_s2,params.RMS ); 
end

% calculate distribution
m       = params.m;
s2      = params.s2;
y       = sqrt(1/(2*pi))./(sqrt(s2)*x).*exp(- ((log(x)-m).^2)/(2*s2));

% plot and write the text
line( 'parent',hAx,'xdata',x,'ydata',y,'linewidth',2,'color',cl );
hTxt    = text( 1.2*mean(x),ylimit(2)*p,txt,'parent',hAx );
ext     = get( hTxt,'Extent' );
line( ext(1) - xlimit/15 ,(ext(2)+ext(4)/2)*[1 1],'linewidth',3,'color',cl );
% ext(3)  = ext(3)+xlimit(2)/15;
% rectangle( 'position',ext );