function [out1,out2] = removerows(in1,in2,in3,in4)
%REMOVEROWS Remove matrix rows with specified indices.
%  
%  Syntax
%
%	[y,ps] = removerows(x,ind)
%	[y,ps] = removerows(fp)
%	y = removerows('apply',x,ps)
%	x = removerows('reverse',y,ps)
%	dx_dy = removerows('dx',x,y,ps)
%	dx_dy = removerows('dx',x,[],ps)
%     name = removerows('name');
%     fp = removerows('pdefaults');
%     names = removerows('pnames');
%     removerows('pcheck',fp);
%
%  Description
%  
%   REMOVEROWS processes matrices by removing rows with the specified indices.
%  
%	REMOVEROWS(X,IND) takes X and an optional parameter,
%	X - NxQ matrix or a 1xTS row cell array of NxQ matrices.
%     IND - Vector of row indices to remove. (Default is [])
%	and returns,
%     Y - Each MxQ matrix, where M==N-length(IND). (optional).
%     PS - Process settings, to allow consistent processing of values.
%
%   REMOVEROWS(X,FP) takes parameters as struct: FP.ind.
%   REMOVEROWS('apply',X,PS) returns Y, given X and settings PS.
%   REMOVEROWS('reverse',Y,PS) returns X, given Y and settings PS.
%   REMOVEROWS('dx',X,Y,PS) returns MxNxQ derivative of Y w/respect to X.
%   REMOVEROWS('dx',X,[],PS)  returns the derivative, less efficiently.
%   REMOVEROWS('name') returns the name of this process method.
%   REMOVEROWS('pdefaults') returns default process parameter structure.
%   REMOVEROWS('pdesc') returns the process parameter descriptions.
%   REMOVEROWS('pcheck',fp) throws an error if any parameter is illegal.
%    
%	Examples
%
%   Here is how to format a matrix so that rows 2 and 4 are removed:
%	
%     x1 = [1 2 4; 1 1 1; 3 2 2; 0 0 0]
%     [y1,ps] = removerows(x1,[2 4])
%
%   Next, we apply the same processing settings to new values.
%
%     x2 = [5 2 3; 1 1 1; 6 7 3; 0 0 0]
%     y2 = removerows('apply',x2,ps)
%
%   Here we reverse the processing of y1 to get x1 again.
%
%     x1_again = removerows('reverse',y1,ps)
%
%  Algorithm
%
%     In the reverse calculation, the unknown values of replaced
%     rows are represented with NaN values.
% 
%  See also MAPMINMAX, FIXUNKNOWNS, MAPSTD, PROCESSPCA, REMOVECONSTANTROWS

% Copyright 1992-2007 The MathWorks, Inc.
% $Revision: 1.1.6.8 $

% Process function boiler plate script
boiler_process

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Name
function n = name
n = 'Remove Rows';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameter Defaults
function fp = param_defaults(values)

num = length(values);
if (num<1), fp.ind =[]; else fp.ind = values{1}; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameter Names
function names = param_names()
names = {'Indices of rows to remove.'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameter Check
function err = param_check(fp)

i = fp.ind;
if ~isa(i,'double') || sum(size(i)> 1)>1 || ~all(isfinite(i)) || ~all(isreal(i)) || any(i<1) || any(i~=floor(i))
  err = 'ind must be a row vector of indices.';
elseif length(unique(i)) ~= length(i)
  err = 'ind must not contain any redundant indices.';
else
  err = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% New Process
function [y,ps] = new_process(x,fp)

[R,Q]=size(x);
if  any(fp.ind > R)
  error('NNET:Process','X has fewer rows than the maximum index to remove.');
end  

% Reduce the transformation matrix appropriately
ps.name = 'removerows';
ps.xrows = R;
ps.yrows = R-length(fp.ind);
ps.remove_ind = fp.ind;
ps.keep_ind = 1:R;
ps.keep_ind(ps.remove_ind) = [];
ps.recreate_ind = zeros(1,R);
ps.recreate_ind(ps.keep_ind) = 1:ps.yrows;

y = apply_process(x,ps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Apply Process
function y = apply_process(x,ps)

y = x(ps.keep_ind,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reverse Process
function x = reverse_process(y,ps)

Q = size(y,2);
x = zeros(ps.xrows,Q);
x(ps.keep_ind,:) = y;
x(ps.remove_ind,:) = NaN; % Don't cares

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Derivative of Y w/respect to X
function dy_dx = derivative(x,y,ps);

Q = size(x,2);
d = zeros(ps.yrows,ps.xrows);
for i=1:length(ps.keep_ind)
  d(i,ps.keep_ind(i)) = 1;
end
dy_dx = d(:,:,ones(1,Q));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Derivative of Y w/respect to X
function dx_dy = reverse_derivative(x,y,ps);

Q = size(x,2);
d = derivative(x(:,1),y(:,1),ps)';
dx_dy = d(:,:,ones(1,Q));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = simulink_params(ps)

p = ...
  { ...
  'inputSize',mat2str(ps.xrows);
  'keep',mat2str(ps.keep_ind);
  };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = simulink_reverse_params(ps)

p = ...
  { ...
  'inputSize',mat2str(ps.xrows);
  'rearrange',mat2str(ps.recreate_ind);
  };
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%