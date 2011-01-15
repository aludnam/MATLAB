%
% This folder contains a collection of "fitting" functions. 
% (Some has demo options - the third section)
% The GENERAL input to the functions should be samples of the distribution.
% 
% for example, if we are to fit a normal distribution ('gaussian') with a mean "u" and varaince "sig"^2
% then the samples will distribute like:  
%      samples = randn(1,10000)*sig + u
% 
%fitting with Least-Squares is done on the histogram of the samples.
% fitting with Maximum likelihood is done directly on the samples.
%
%
% Contents of this folder
% =======================
% 1) Maximum likelihood estimators
% 2) Least squares estimators
% 3) EM algorithm for estimation of multivariant gaussian distribution (mixed gaussians)
% 4) added folders: Create - which create samples for the EM algorithm test
%                   Plot   - used to plot each of the distributions (parametric plot)
%
%
%
%
%
% Maximum likelihood estimators
% =============================
% fit_ML_maxwell   - fit maxwellian distribution
% fit_ML_rayleigh  - fit rayleigh distribution
%                      (which is for example: sqrt(abs(randn)^2+abs(randn)^2))
% fit_ML_laplace   - fit laplace distribution
% fit_ML_log_normal- fit log-normal distribution
% fit_ML_normal    - fit normal (gaussian) distribution
%
% NOTE: all estimators are efficient estimators. for this reason, the distribution
%       might be written in a different way, for example, the "Rayleigh" distribution
%       is given with a parameter "s" and not "s^2".
%
%
% least squares estimators
% =========================
% fit_maxwell_pdf  - fits a given curve of a maxwellian distribution
% fit_rayleigh_pdf - fits a given curve of a rayleigh distribution
%
% NOTE: these fit function are used on a histogram output which is like a sampled 
%       distribution function. the given curve MUST be normalized, since the estimator
%       is trying to fit a normalized distribution function.
%
%
%
%
% Multivariant Gaussian distribution
% ==================================
% for demo of 1D mixed-gaussian fitting, run:  fit_mix_gaussian
% for demo of 2D mixed-gaussian fitting, run:  fit_mix_2d_gaussian
%
% these routines fit and plot the results of the parameters of: 
% random distribution of random amount of gaussians with random parameters 
%
%
