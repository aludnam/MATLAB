function x_fit = select_fitrange(x, fitrange)
% SELECT_FITRABGE selects range for fitting
% x_fit = select_fitrange(x, fitrange)
% x - valueas on the x-axis
% fitrange = [low_limit, high_limit] - range (in units of the x-axis) to select 

x_fit_bin = and(x > fitrange(1), x < fitrange(2));
x_fit = x .* x_fit_bin;