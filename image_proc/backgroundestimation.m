function [out_nobg, peval]=backgroundestimation(dpix, peval)
% [out_nobg, peval]=backgroundestimation(dpix, peval)

dpix_dip = dip_image(dpix);
% Assigning some empirical default values
if ~isfield(peval, 'bg_clip'); peval.bg_clip = 'no'; end
if ~isfield(peval, 'bg_fs_var'); peval.bg_fs_var = 5; end
if ~isfield(peval, 'bg_perc'); peval.bg_perc = 20; end
if ~isfield(peval, 'bg_ob_dist'); peval.bg_ob_dist = 8; end

[out_nobg, peval.bg, bg_im]=backgroundoffset(dpix_dip, peval.bg_clip, peval.bg_fs_var, peval.bg_perc, peval.bg_ob_dist);