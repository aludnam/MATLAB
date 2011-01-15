function plot_maxwell( x,params,hAx,plot_num,fontsize )
% plot the maxwell distribution with parameter "a"
% 
% the distribution is given by:
%
%        p(r) = sqrt(2/pi)*(a^(-3/2))*(r^2)*exp(-(r^2)/(2*a))
%
% format:   plot_maxwell( x,params,hAx,plot_num,fontsize )
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
% example:  plot_maxwell( x,fit_ML_maxwell(data),hAx,3 )
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
fnc_txt = '\surd\it{2/\pi}  \bf{\cdot a}^{-3/2}\bf{x}^2\bf{e}^{-x^2/2a}\rm';
if isfield( params,'VAR' )
    txt     = sprintf( '\\fontsize{%d}\\bfMaxwell PDF\\rm with %s:  %s\na               = %g\nVAR(a)      = %1.3g\nRMS err(a) = %1.3g\n',...
        fontsize,params.type,fnc_txt,params.a,params.VAR,params.RMS );
else
    txt     = sprintf( '\\fontsize{%d}\\bfMaxwell PDF \\rmwith %s:  %s\na               = %g\nCRB(a)      = %1.3g\nRMS err(a) = %1.3g\n',...
        fontsize,params.type,fnc_txt,params.a,params.CRB,params.RMS ); 
end

% calculate distribution
a       = params.a;
x2      = x.^2;
y       = sqrt(2/pi)*x2.*a^(-1.5).*exp(-x2/(2*a));

% plot and write the text
line( 'parent',hAx,'xdata',x,'ydata',y,'linewidth',2,'color',cl );
hTxt    = text( 1.2*mean(x),ylimit(2)*p,txt,'parent',hAx );
ext     = get( hTxt,'Extent' );
line( ext(1) - xlimit/15 ,(ext(2)+ext(4)/2)*[1 1],'linewidth',3,'color',cl );
% ext(3)  = ext(3)+xlimit(2)/15;
% rectangle( 'position',ext );