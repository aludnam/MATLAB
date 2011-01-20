function [varargout, peval] = separcomp(dpix, peval, winit_pix, hinit, verbose)
% separcomp(dpix, p_script, savethis, winit, hinit, verbose)
% separate components
% V ~ WH
% V -> N_pix x N_t
% W -> N_pix x N_comp - xth pixel of the ith components
% H -> N_copm x N_t - contribution of the i-th component in the time t
tstart = tic;
if ~exist('verbose','var')
    verbose = 1;
end

if ~isfield(peval, 'w_fixvec') peval.w_fixvec=[]; end
if ~isfield(peval, 'h_fixvec') peval.h_fixvec=[]; end

% Fixing background components:
if ~isfield(peval, 'fixbackground') 
    peval.fixbackground = 'wh';
end

% Fixing w component of the background
if sum(peval.fixbackground == 'w')>0 
    peval.w_fixvec = [peval.w_fixvec, peval.ncomp + 1];
end

% Fixing h component of the background
if sum(peval.fixbackground == 'h')>0 
    peval.h_fixvec = [peval.h_fixvec, peval.ncomp + 1];
end

[peval.nx, peval.ny, peval.nt] = size(dpix);
peval.numpix = peval.nx*peval.ny;
peval.meandata = mean(dpix(:));

dvec = reshape(dpix,peval.numpix, peval.nt);

% Background subtraction with empirical values:
if ~isfield(peval, 'bg')
    [out_nobg, peval]=backgroundestimation(dpix, peval);
end

% Estimation of the number of componenets
if ~isfield(peval, 'ncomp')
    peval.ncomp=estimate_ncpomp(dvec);
end

% Initialization of w and h
[winit, hinit] = initwh(winit_pix, hinit, peval, verbose);

% This is to fix zeros in the data as the divergence is not defined for that...
dvec=max(dvec, eps);

eval (['[w,h,peval, dtrace, htrace]=' peval.method '(dvec,winit,hinit,peval,verbose);']);
peval.ddiv_end = ddivergence(dvec, w*h);
peval.gradddiv_end = dtrace(end-1)-dtrace(end);
%better with feval....
% varargout = struct('w',w,'h',h, 'dtrace', dtrace, 'htrace', htrace);
varargout = struct('w',w,'h',h);
peval.elapsedtime_separcomp = toc(tstart);