function plot_rayleigh( x,params,hAx,plot_num,fontsize )
% plot the rayleigh distribution with parameter "a"
% 
% the distribution is given by:
%
%        p(r)=r*exp(-r^2/(2*s))/s
%
% format:   plot_rayleigh( x,params,hAx,plot_num,fontsize )
%
% input:    x         - X axis, for the plot
%           params    - the distribution parameter, RMS error and VAR or CRB
%           hAx       - where to plot the distribution curve
%           plot_num  - since the curve is added with a text to the axes,
%                       this parameter specifies where the text should be displayed
%                       and what color to choose for the curve
%           fontsize  - size of the font of the text, default 9
%
% example:  plot_rayleigh( x,fit_ML_rayleigh(data),hAx,3 )
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
fnc_txt = '(\bf{x/s}\rm)\cdot\bfe\rm^{-\bfx^2/2\bfs}\rm';
if isfield( params,'VAR' )
    txt     = sprintf( '\\fontsize{%d}\\bfRayleigh PDF\\rm with %s:  %s\ns               = %g\nVAR(s)      = %1.3g\nRMS err(s) = %1.3g\n',...
        fontsize,params.type,fnc_txt,params.s,params.VAR,params.RMS );
else
    txt     = sprintf( '\\fontsize{%d}\\bfRayleigh PDF\\rm with %s:  %s\ns               = %g\nCRB(s)      = %1.3g\nRMS err(s) = %1.3g\n',...
        fontsize,params.type,fnc_txt,params.s,params.CRB,params.RMS );
end

% calculate distribution
s       = params.s;
y       = x.*exp(-x.^2/(2*s))/s;

% plot and write the text
line( 'parent',hAx,'xdata',x,'ydata',y,'linewidth',2,'color',cl );
hTxt    = text( 1.2*mean(x),ylimit(2)*p,txt,'parent',hAx );
ext     = get( hTxt,'Extent' );
line( ext(1) - xlimit/15 ,(ext(2)+ext(4)/2)*[1 1],'linewidth',3,'color',cl );
% ext(3)  = ext(3)+xlimit(2)/15;
% rectangle( 'position',ext );