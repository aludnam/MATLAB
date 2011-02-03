function [out_nobg, peval]=backgroundestimation(dpix, peval, bgestimate)
% [out_nobg, peval]=backgroundestimation(dpix, peval, clip, bgestimate)
% bgestimat: (optional) if supplied then used as a background
% clip: (default 'no') clip negative values

if ~isfield(peval, 'bg_clip'); peval.bg_clip = 'no'; end
if exist('bgestimate', 'var')
    out_nobg = dpix-bgestimate;
    peval.bg = bgestimate;
    msg = 'Background supplied.';
else
    
    dpix_dip = dip_image(dpix);
    % Assigning some empirical default values
    if ~isfield(peval, 'bg_fs_var'); peval.bg_fs_var = 5; end
    if ~isfield(peval, 'bg_perc'); peval.bg_perc = 20; end
    if ~isfield(peval, 'bg_ob_dist'); peval.bg_ob_dist = 8; end
    
    [out_nobg, peval.bg, bg_im]=backgroundoffset(dpix_dip, peval.bg_clip, peval.bg_fs_var, peval.bg_perc, peval.bg_ob_dist);
    msg = 'Background estimated using ''backgroundoffset'' function.';
end
mfprintf(peval.fid, [msg '(background = %g)\n'], peval.bg)
if strcmpi (peval.bg_clip, 'yes') 
    out_nobg = max(out_nobg, eps);
    mfprintf(peval.fid, 'Negative and zero values clipped when subtracting background.\n');
end
    
