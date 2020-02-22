%FINDSPECIES  Find indices for species in a model or compartment
%
% findspecies(modelObj, namelist) returns a list of indices for the
% species given in namelist.  The names should be specified as a
% cell-array.

% RMM, 7 Sep 2012

function indexlist = findspecies(modelObj, namelist,varargin)

listOfSpecies = get(modelObj.species, 'name');
if ischar(listOfSpecies)
    listOfSpecies = {listOfSpecies};
end

if ischar(namelist)
    namelist = {namelist};
end


if nargin > 2
    if strcmpi(varargin{1},'withInComplex')
        indexlist = cellfun(@(x) findStringInAList(listOfSpecies,x,'withinText'),namelist,'UniformOutput',false);
        indexlist = cell2mat(indexlist);
    end
else
    indexlist = cellfun(@(x) findStringInAList(listOfSpecies,x),namelist);
end


% Local variables:
% mode: matlab
% End:
