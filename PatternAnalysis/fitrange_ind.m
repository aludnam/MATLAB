function [ind_low, ind_high] = fitrange_ind(x, fitrange)
% SELECT_FITRABGE gives indices in x corresponding to values in fitrange
% [ind_low, ind_high] = fitrange_ind(x, fitrange)
% x - valueas on the x-axis
% fitrange = [low_limit, high_limit] - range (in units of the x-axis) to select 

x_ind = find(and(x > fitrange(1), x < fitrange(2)));
ind_low = x_ind(1);
ind_high = x_ind(end);