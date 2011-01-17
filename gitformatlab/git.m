function git(varargin)
%   A thin MATLAB wrapper for the git source control system
%   The short instructions are:
%       Use this exactly as you would use the OS command-line Git.
% 
%   The long instructions are:
%       This is not meant to be a comprehensive guide to near-omnipotent 
%       Git source control. That is here:
%           http://git-scm.com/documentation
% 
%       Common MATLAB workflow: 
%    
%       % Creates initial repository tracking all files under some root
%       % folder
%       >> cd ~/
%       >> git init               
%     
%       % Shows changes made to all files in repo (none so far)
%       >> git status             
%    
%       % Create a new file and add some code
%       >> edit foo.m       
% 
%       % Check repo status, after new file created
%       >> git status
%      
%       % Add foo.m to your repo
%       >> git add foo.m   
%     
%       % Commit your changes to a new branch, with comments
%       >> git commit -m 'Created new file, foo.m'
% 
%       % Other useful commands (replace ellipses with appropriate args)
%       >> git checkout ...       % To return to an earlier node in branch
%       >> git branch ...         % To create or move to another branch
%       >> git merge ...          % When merge conflicts arise
%       >> git diff ...           % See line-by-line changes 
%    
%   Useful resources:
%       1. GitX: A visual interface for Git on the OS X client
%       2. Github.com: Remote hosting for Git repos
%       3. Git on Wikipedia: Further reading 
% 
% v0.1,     27 October, 2010 -- Initial support for OS X & Linux, untested
%                               on PCs, but expected to work
%
% Author: Manu Raghavan
% Free for commercial use

% if ismac || isunix
%     [status] = system('git');
%     assert(status==1, ...
%         strcat('git is not installed\n',...
%                'Download it at http://git-scm.com/download'));
% elseif ispc
%     warning('Make sure git is installed');
% end
    
arguments = parse(varargin{:});  
[smaz,result] = system(['git ',arguments]);
disp(result)
end

function space_delimited_list = parse(varargin)
    space_delimited_list = cell2mat(...
                cellfun(@(s)([s,' ']),varargin,'UniformOutput',false));
end