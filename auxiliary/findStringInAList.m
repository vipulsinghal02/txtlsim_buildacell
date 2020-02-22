% finds a string match (based on strcmp) in a string list and gives back
% its position(s) in the list.
% input: list of string, string of interest
% output: all the exact matches' index

function varargout = findStringInAList (list,string,varargin)

% selecting operation mode
if nargin > 2
   if strcmpi(varargin{1},'withinText');
    resList = strfind(list,string);
    indexList = find(~cellfun('isempty',resList) >0);
   else 
     error('not valied option: %s',varargin{1});
   end
else
    
    if ~isempty(list) && ~isempty(string)
        if ischar(list)
            indexList = strcmp(list,string);
        else
            indexes = cellfun(@(x) strcmp(x,string),list);
            ind = find(indexes > 0);
            if ~isempty(ind)
                indexList = ind;
            else
                indexList = 0;
            end
        end
    else
        indexList = 0;
    end
end

% building output
varargout{1} = indexList;
if nargout >1
    varargout{2} = string;
end