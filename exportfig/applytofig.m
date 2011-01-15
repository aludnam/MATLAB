function old = applytofig(h,varargin)
%APPLYTOFIG  Apply EXPORTFIG options to figure.
%   OLD = APPLYTOFIG(H) applies the default EXPORTFIG options to
%   the figure with handle H and returns the state-difference
%   structure in OLD.
%   OLD = APPLYTOFIG(H, OPTIONS) applies OPTIONS to H as described
%   in EXPORTFIG.
%   OLD = APPLYTOFIG(...,PARAM1,VAL1,PARAM2,VAL2,...) applies the
%   specified parameter-value pairs to H as described in EXPORTFIG.
%   
%   Use RESTOREFIG(H,OLD) to restore the figure to its state before
%   applyfig was called.
%
%   See also EXPORTFIG, PREVIEWFIG, RESTOREFIG.

%  Copyright 2000-2009 The MathWorks, Inc

temp = tempname;
old = exportfig(h,temp, varargin{:},'applystyle',1);
