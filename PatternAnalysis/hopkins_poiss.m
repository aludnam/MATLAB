function [Hpoiss] = hopkins_poiss(xH, m)
% HOPKINS_POISS plots hopkins statistics for Poisson process
% [Hpoiss] = hopkins_poiss(xH, m)
% xH - position on x-axis where Hopkins is computed
% m - parameter of the hopkins statistics

Hpoiss = gamma(2*m)/(gamma(m)^2) * xH.^(m-1) .* (1-xH).^(m-1);