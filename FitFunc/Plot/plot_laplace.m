function plot_laplace( x,params,hAx,plot_num,fontsize )
% plot the laplace distribution with parameter "u" and "b"
% 
% the distribution is given by:
%
%        p(x) = 1/(2*b)*exp(-abs(x-u)/b)
%
% format:   plot_laplace( x,params,hAx,plot_num,fontsize )
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
% example:  plot_laplace( x,fit_ML_laplace(data),hAx,3 )
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
fnc_txt = '\it{1/2\bfb}\rm \bf{\cdot e}^{-\mid\bf{x-\mu}\mid/\bf{b}}\rm';
if isfield( params,'VAR' )
    txt     = sprintf( '\\fontsize{%d}\\bfLaplace PDF\\rm with %s:  %s\n\\mu = %g    b = %g\nVAR(b) = %1.3g\nRMS err = %1.3g\n',...
        fontsize,params.type,fnc_txt,params.u,params.b,params.VAR_b,params.RMS );
else
    txt     = sprintf( '\\fontsize{%d}\\bfLaplace PDF\\rm with %s:  %s\n\\mu = %g    b = %g\nCRB(b) = %1.3g\nRMS err = %1.3g\n',...
        fontsize,params.type,fnc_txt,params.u,params.b,params.CRB_b,params.RMS ); 
end

% calculate distribution
u       = params.u;
b       = params.b;
y       = 1/(2*b)*exp(-abs(x-u)/b);

% plot and write the text
line( 'parent',hAx,'xdata',x,'ydata',y,'linewidth',2,'color',cl );
hTxt    = text( 1.2*mean(x),ylimit(2)*p,txt,'parent',hAx );
ext     = get( hTxt,'Extent' );
line( ext(1) - xlimit/15 ,(ext(2)+ext(4)/2)*[1 1],'linewidth',3,'color',cl );
% ext(3)  = ext(3)+xlimit(2)/15;
% rectangle( 'position',ext );