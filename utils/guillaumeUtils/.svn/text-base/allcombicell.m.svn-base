function varargout = allcombi(varargin)
% ALLCOMBICELL - All the cell combinations
%
% [C1,C2] = allcombi(A1,A2) 
%   Ai is a cell array
%   Ci is a cell array
%
% example:
% --------
% [C1,C2]=allcombicell({'err','ds'},{'jlkmq','jkl'})
% returns
% C1 = 
%     'err'    'err'    'ds'    'ds'
% C2 = 
%     'jlkmq'    'jkl'    'jlkmq'    'jkl'
%
%
% 
% See also ALLCOMBI, MESHGRID

% To have a result in a cell array of cells:
% C=cell(1,2);[C{:}]=allcombicell({'err','ds'},{'jlkmq','jkl'})

ic = find(cellapply(varargin,'~iscell'));%arguments that are not cells
varargin(ic) = num2cell(varargin(ic)); %make them cells

L = cellfun('length',varargin);
if any(L==0)
    varargout = {{}};
    return
end


ac = allcombi(L);

if nargout<=1
    nc = size(ac,1);    
    varargout = cell(nc,nargin);
    for i=1:nargin
        [varargout{:,i}] = deal(varargin{i}{ac(:,i)'});
    end
    varargout={num2cell(varargout,2)};
else
    varargout = cell(1,nargin);
    for i=1:nargin
        varargout{i} = varargin{i}(ac(:,i)');
    end

    if nargout<length(varargout)
        error('not enough arguments');
    end
end
