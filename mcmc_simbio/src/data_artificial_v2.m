function [di] = data_artificial_v2(mobj, tv, measuredSpecies, ...
    dosedNames, dosedVals, activeNames, activeValues, varargin)
% [di] = data_artificial_v2(mobj, tv, measuredSpecies, ...
%     dosedNames, dosedVals, activeNames, activeValues, varargin)
% use simbiology model object to generate artificial data.
%
% There are two ways to define the inputs. In the SCALAR MODE, we have:
%
% mobj: simbiology model object
%
% tv: a vector of timepoints which the data in the output data array will
% correspond to. These are the poitns at which the model will be simulated to
% compute the data values.
%
% measuredSpecies: this is a cell array of cell arrays of strings. The
% species in the inner cell array get summed to create the measured
% species. The order of the outer cell array defines the ordering of the
% second dimension of the dataArray property of the data_info struct that
% is output.
%
% dosedNames: Cell array of species that get dosed.
%
% dosedVals: matrix of size #number of dosed species x number of dose
% combinations.
%
% activeNames: cell array of parameters and species names to set in the
% model.
%
% activeValues: vector of corresponding values. Note that the parameter
% values are NOT log transformed, ie, they all lie in the nonnegative orthant.
%
% In the CELL MODE, the main difference is that the data_info output struct
% an now be a non singleton array of length nDatasets. The activeValues
% input is now a cell array of numerical vectors. It has dimensions
% nDatasets x 1 or 1 x nDatasets.
% 
% mobj: A 1x1 cell containing a simbiology model object, or a row or column
% cell array of length nDatasets containing model objects.
%
% tv: A 1x1 cell containing a vector of timepoints, or a row or column
% cell array of length nDatasets containing vectors of timepoints.  See
% scalar version documentation above for more info.
%
% measuredSpecies: A 1x1 cell containing a cell array of cell arrays of
% strings, or a row or column cell array of length nDatasets containing
% cell array of cell arrays of strings. See scalar version documentation
% above for more info.
%
% dosedNames: A 1x1 cell containing a Cell array of species that get dosed,
% or a row or column cell array of length nDatasets containing
% cell arrays of species that get dosed. See scalar version documentation
% above for more info.
%
% dosedVals: A 1x1 cell containing a matrix of size
% #number of dosed species x number of dose combinations,
% or a row or column cell array of length nDatasets containing
% matrix of size #number of dosed species x number of dose combinations.
% See scalar version documentation above for more info.
%
% activeNames: A 1x1 cell containing a cell array of parameters and species
% names to set, or a row or column cell array of length nDatasets containing
% cell arrays of parameters and species names to set. See scalar version
% documentation above for more info.
%
% activeValues: A 1x1 cell containing a vector of corresponding values,
% or a row or column cell array of length nDatasets containing
% vectors of corresponding values. See scalar version documentation
% above for more info.
%
%
% OPTIONAL NAME VALUE PAIR ARGUMENTS - used to populate the data_info struct.
% If cell mode is active, then all of these are correspindingly encapsulated
% in a 1x1 cell or a row or column cell array of length nDatasets.
%
% 'dataInfo': Human readable description of the data. Should be specified using
% the format that makes it printable with the function fprintf. (so newline
% characters \n at the ver least, for example.). If not specified, its value is
% simply 'Artificial Data'.
%
% 'timeUnits': String specifying the units of the time axis. Default is 'seconds'
%
% 'dataUnits': A cell array of strings specifying the units of the measured
% species. If not specified, the default is 'nM'.
%
% 'doseUnits': A cell array of strings specifying the units of the
% dosed species. If not specified, the default is 'nM'.
%
% 'noise', VALUE; where VALUE is a vector of standard deviations of the gaussian
% noise added to the data, each element corresponding to one measured
% species.
%
% 'replicates', VALUE, where VALUE is a positive interger of the number
% of replicates.
%
%
% OUTPUTS: This function returns a data_info struct with fields
% ------------------------------------------
%
% 'dataInfo': A human readable text description of the data. If not specified
% as a name value pair input argument, the string 'Artificial Data' is used.
%
% 'timeVector': vector of timepoints, same as tv, a required positional input.
%
% 'timeUnits': units of the time vector. Allowed options are:
% 'seconds', 'minutes', 'hours', 'days', 'weeks'. If no units are specified
% as a name value input, then the units are specified as 'seconds'.
%
% 'dataArray': An array contianing the raw data that is generated by simulating
% the data according to the mcmc_info struct. Typically has dimensions
% corresponding to timepoints x measured outputs x replicates x doses.
%
% 'measuredNames': A 1 x number of measured species cell array of the strings
% specifying which species are dosed. These are not strings corresponding to
%  the species in the model. Takes from the corresponding name value pair input
% argument. If not specified, the values are taken from the mcmc_info struct.
%
% 'dataUnits': A 1 x number measured species cell array of units corresponding to
% the raw data in the dataArray. If no units are specified, then nM are used.
%
% 'dimensionLabels': a 1 by length(size(data_info.dataArray)) cell array of
% labels for the dimensions of the dataArray.
%
% 'dosedNames': A 1 x number of dosed species cell array of the strings specifying
% which species are dosed. These are not strings corresponding to the species
% in the model. See mcmc_info constructor functions for that.
%
% 'dosedVals': A matrix of dose values of size
% # of dosed species by # of dose combinations
%
% 'doseUnits': A 1 x number of dosed species cell array of strings specifying the
% units of the dosed species. If no units are
%

% Copyright (c) 2018, Vipul Singhal, Caltech
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.


% check the sizes of things:
if iscell(activeValues)
    nDatasets = length(activeValues);
    cellMode = true;
else
    nDatasets = 1;
    cellMode = false;
end

if nDatasets == 1
    assert(~iscell(mobj) && length(mobj) == 1);
    assert(~iscell(tv) );
    assert(~iscell(dosedVals));
end

if cellMode
    
    p = inputParser;
    
    addParameter(p, 'replicates', {1});  %
    addParameter(p, 'noise', 'default');  % the default is 0
    addParameter(p, 'timeUnits', {'seconds'});
    addParameter(p, 'dataUnits', {'nM'});
    addParameter(p, 'doseUnits', {'nM'});
    addParameter(p, 'dataInfo', {'Artificial Data'});
    addParameter(p, 'dimensionLabels', {{'time points', 'measured species',...
        'replicates', 'doses'}});
    parse(p, varargin{:});
    p = p.Results;
    for i = 1:nDatasets
        currmeasuredSpecies = cellcontents(measuredSpecies, i);
        
        if strcmp(p.noise, 'default')
            noisevec = zeros(length(currmeasuredSpecies), 1);
        else
            noisevec = cellcontents(p.noise, i);
        end
        nReplicates = cellcontents(p.replicates, i);
        timeUnits = cellcontents(p.timeUnits, i);
        dataUnits = cellcontents(p.dataUnits, i);
        doseUnits = cellcontents(p.doseUnits, i);
        dataInfo = cellcontents(p.dataInfo, i);
        dimensionLabels = cellcontents(p.dimensionLabels, i);
        currmobj = cellcontents(mobj, i);
        currtv = cellcontents(tv, i);
        currdosedNames = cellcontents(dosedNames, i);
        currdosedVals = cellcontents(dosedVals, i);
        curractiveNames = cellcontents(activeNames, i);
        if ~length(activeValues)>1
            error('what on earth?')
        end
        curractiveValues = activeValues{i};
        
        da = computeArtificialData(currmobj, currtv, noisevec, ...
            currmeasuredSpecies, nReplicates, currdosedVals, currdosedNames, ...
            curractiveNames, curractiveValues);
        
        if i ==1
            % make the data_info struct
            di = struct('dataInfo', {dataInfo}, ...
                'timeVector', {currtv}, ...
                'timeUnits', {timeUnits},...
                'dataArray', {da},...
                'measuredNames', {currmeasuredSpecies},...
                'dataUnits', {dataUnits},...
                'dimensionLabels', {dimensionLabels}, ...
                'dosedNames', {currdosedNames},...
                'dosedVals', {currdosedVals}, ...
                'doseUnits', {doseUnits});
        else
            currdi = struct('dataInfo', {dataInfo}, ...
                'timeVector', {currtv}, ...
                'timeUnits', {timeUnits},...
                'dataArray', {da},...
                'measuredNames', {currmeasuredSpecies},...
                'dataUnits', {dataUnits},...
                'dimensionLabels', {dimensionLabels}, ...
                'dosedNames', {currdosedNames},...
                'dosedVals', {currdosedVals}, ...
                'doseUnits', {doseUnits});
            di = [di; currdi];
        end
    end
else
    noisedefaults = zeros(length(measuredSpecies), 1);
    dimensionLabels = {'time points', 'measured species', 'replicates', 'doses'};
    p = inputParser;
    addParameter(p, 'replicates', 1, @isnumeric);  %
    addParameter(p, 'noise', noisedefaults, @isnumeric);  % the default is 0
    addParameter(p, 'timeUnits', 'seconds', @ischar);
    addParameter(p, 'dataUnits', 'nM');
    addParameter(p, 'doseUnits', 'nM');
    addParameter(p, 'dataInfo', 'Artificial Data');
    addParameter(p, 'dimensionLabels', dimensionLabels);
    parse(p, varargin{:});
    p = p.Results;
    noisevec = p.noise;
    nReplicates = p.replicates;
    
    da = computeArtificialData(mobj, tv, noisevec, measuredSpecies, ...
        nReplicates, dosedVals, dosedNames, activeNames, activeValues);
    
    % make the data_info struct
    di = struct('dataInfo', {p.dataInfo}, ...
        'timeVector', {tv}, ...
        'timeUnits', {p.timeUnits},...
        'dataArray', {da},...
        'measuredNames', {measuredSpecies},...
        'dataUnits', {p.dataUnits},...
        'dimensionLabels', {p.dimensionLabels}, ...
        'dosedNames', {dosedNames},...
        'dosedVals', {dosedVals}, ...
        'doseUnits', {p.doseUnits});
    
end
end

function da = computeArtificialData(mobj, tv, noisevec, measuredSpecies, ...
    nReplicates, dosedVals, dosedNames, activeNames, activeValues)

configsetObj = getconfigset(mobj, 'active');
set(configsetObj, 'StopTime', tv(end));


% initialize the data array
da = zeros(length(tv), length(measuredSpecies),...
    nReplicates, size(dosedVals, 2));

% set parameters and species initial concentrations in the model
for i = 1:length(activeNames)
    p1 = sbioselect(mobj.parameters, 'Name',  activeNames{i});
    if ~isempty(p1)
        set(p1, 'Value', activeValues(i));
    end
end

for i = 1:length(activeNames)
    s1 = sbioselect(mobj.species, 'Name', activeNames{i});
    if ~isempty(s1)
        set(s1, 'InitialAmount', activeValues(i))
    end
end

% set dose values, simulate model, and populate output data array.
for dID = 1:size(dosedVals, 2)
    % set the dose value using the mcmc_info struct
    for i = 1:length(dosedNames)
        s1 = sbioselect(mobj.species, 'Name', dosedNames{i});
        if ~isempty(s1)
            set(s1, 'InitialAmount', dosedVals(i, dID))
        end
    end
    
    % simulate the model.
    sd = sbiosimulate(mobj);
    sd = resample(sd, tv);
    for msID = 1:length(measuredSpecies)
        currmeasuredSpecie = measuredSpecies{msID};
        spSD = selectbyname(sd, currmeasuredSpecie);
        summed_trajectories = sum(spSD.Data, 2);
        %  add noise if needed.
        for rID = 1:nReplicates
            da(:, msID, rID, dID) = ...
                summed_trajectories + ...
                noisevec(msID)*randn(length(tv), 1);
        end
    end
end


end

function currcellcontents = cellcontents(cellarray, count)
assert(iscell(cellarray));

if length(cellarray)>1
    currcellcontents = cellarray{count};
else
    currcellcontents = cellarray{1};
end

end



