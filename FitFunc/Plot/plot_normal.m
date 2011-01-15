function plot_normal( x,params,hAx,plot_num,fontsize )
% plot the normal distribution with parameter "u" and "sig2"
% 
% the distribution is given by:
%
%        p(r) = sqrt(1/2/pi/sig^2)*exp(-((r-u)^2)/(2*sig^2))
%
% format:   plot_normal( x,params,hAx,plot_num,fontsize )
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
% example:  plot_normal( x,fit_ML_normal(data),hAx,3 )
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

% calculate distribution
u       = params.u;
sig2    = params.sig2;
y       = sqrt(1/2/pi/sig2)*exp(-(x-u).^2/(2*sig2));

% plot and get axis limits
line( 'parent',hAx,'xdata',x,'ydata',y,'linewidth',2,'color',cl );
ylimit  = ylim(hAx);
xlimit  = xlim(hAx);

% decide where to put the text, and the content of the text
p       = plot_num*0.15 + 0.3;
fnc_txt = '\surd\it{1/2\pi \sigma^2}\bf{\cdot e}^{-(x-\mu)^2/2\bf{\sigma^2}}\rm';
if isfield( params,'VAR' )
    txt     = sprintf( '\\fontsize{%d}\\bfNormal PDF\\rm with %s:  %s\n\\mu = %g    \\sigma^2 = %g\nVAR(\\mu) = %1.3g  VAR(\\sigma^2) = %1.3g\nRMS err = %1.3g\n',...
        fontsize,params.type,fnc_txt,params.u,params.sig2,params.VAR_u,params.VAR_sig2,params.RMS );
else
    txt     = sprintf( '\\fontsize{%d}\\bfNormal PDF\\rm with %s:  %s\n\\mu = %g    \\sigma^2 = %g\nCRB(\\mu) = %1.3g  CRB(\\sigma^2) = %1.3g\nRMS err = %1.3g\n',...
        fontsize,params.type,fnc_txt,params.u,params.sig2,params.CRB_u,params.CRB_sig2,params.RMS ); 
end

% write the text
hTxt    = text( 1.2*mean(x),ylimit(2)*p,txt,'parent',hAx );
ext     = get( hTxt,'Extent' );
line( ext(1) - [0 max(xlimit/15)] ,(ext(2)+ext(4)/2)*[1 1],'linewidth',3,'color',cl );
% ext(3)  = ext(3)+xlimit(2)/15;
% rectangle( 'position',ext );