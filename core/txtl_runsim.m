%
%
% Written by Zoltan A Tuza and Vipul Singhal, Sep 2012
%
% Copyright (c) 2012 by California Institute of Technology
% All rights reserved.
%
%
%
%%%%%% DEFAULT MODE %%%%%%
% [t_ode, x_ode, modelObj, simData] = txtl_runsim(modelObj, configsetObj, time_vector, data, simData)
%
% Input combinations:
% modelObj (this runs a parameter estimation mode, no simulation)
% modelObj, configsetObj
% modelObj, configsetObj, simData
% modelObj, configsetObj, time_vector, data
%
%
% Output combinations:
% simData
% t_ode, x_ode
% t_ode, x_ode, modelObj
% t_ode, x_ode, modelObj, simData
%
%%%%%% EVENTS MODE %%%%%%
% [t_ode, x_ode, modelObj, simData] = txtl_runsim(modelObj, configsetObj, eventTriggers, eventFcns, time_vector, data, simData, 'events')
%
% Input combinations:
% modelObj, configsetObj, eventTriggers, eventFcns, 'events'
% modelObj, configsetObj, eventTriggers, eventFcns, time_vector, data, 'events'
% modelObj, configsetObj, eventTriggers, eventFcns, time_vector, data, simData, 'events'
%
%
% when using this, it is Necessary to use the model object returned by this
% function in subsequent calls / plotting. This is because the first call to
% this function sets the reactions in the modelObject.
%

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%   1. Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%
%   2. Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%      documentation and/or other materials provided with the distribution.
%
%   3. The name of the author may not be used to endorse or promote products
%      derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
% INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
% STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
% IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.



function [varargout] = txtl_runsim(varargin)
%%


modelObj = varargin{1};
if nargin >1
    if isnumeric(varargin{2})
        configsetObj = getconfigset(modelObj, 'active');
        simulationTime = varargin{2};
        set(configsetObj, 'SolverType', 'ode15s');
        set(configsetObj, 'StopTime', simulationTime);
        
        
    else
        configsetObj = varargin{2};
    end
end
if ~verLessThan('Matlab','8')
    configsetObj = getconfigset(modelObj, 'active');
    Sopt = get(configsetObj,'SolverOptions');
    set(Sopt,'AbsoluteToleranceScaling',false);
end


switch nargin
    case 1
        % parameter estimation mode, runsim just assemble the
        % reactions, but won't run it
        configsetObj = [];
        time_vector = [];
        data =  [];
    case 2
        time_vector = [];
        data =  [];
    case 3
        simData = varargin{3};
        [time_vector, prelimData, metaData] = getdata(simData, 'nummetadata');
        species_ind = cellfun(@(x) strcmp(x.Type,'species'), metaData);
        data = prelimData(:,species_ind);
    case 4
        time_vector = varargin{3};
        data = varargin{4};
        
        %         case 7
        %
        %             eventTriggers = varargin{3};
        %             eventFcns = varargin{4};
        %             time_vector = varargin{5};
        %             data = varargin{6};
        %         case 8
        %
        %             eventTriggers = varargin{3};
        %             eventFcns = varargin{4};
        %             time_vector = varargin{5};
        %             data = varargin{6};
        %             simData = varargin{7};
    otherwise
        error('txtl_runsim should be called with 1 to 5 arguments.');
end

m = get(modelObj, 'UserData');

% check what proteins are present, but no corresponding DNA. this will mean
% that the protein must have been added. So now set up the reactions for
% that protein. and you are done! So basically, you are comparing two
% lists.

add_dna_mode = struct('add_dna_driver', {'Setup Reactions'});
txtl_setup_new_protein_added(modelObj, add_dna_mode);

%%% handling events

%     evt = cell(1, length(eventTriggers));
%     for i = 1:length(eventTriggers)
%         evt{i} = addevent(modelObj, eventTriggers{i}(1:end), eventFcns{i}(1:end));
%     end



%% FIRST RUN
% if m.FirstRun
%     for i = 1:length(m.DNAinfo) % should we not set up reactions again if they have already been set up
%         %(we reset the firstRun flag if more dna is added, and then the reactions for ALL the DNA are resetup.)
%         txtl_add_dna(modelObj, m.DNAinfo{i}{1}, m.DNAinfo{i}{2}, ...
%             m.DNAinfo{i}{3}, m.DNAinfo{i}{4}, m.DNAinfo{i}{5}, 'Setup Reactions');
%     end
%     m.FirstRun = False;
% end


for i = 1:length(m.DNAinfo) % should we not set up reactions again if they have already been set up
    %(we reset the firstRun flag if more dna is added, and then the reactions for ALL the DNA are resetup.)
    if strcmp(m.DNAinfo{i}{6}, 'rxns_not_set_up')
        txtl_add_dna(modelObj, m.DNAinfo{i}{1}, m.DNAinfo{i}{2}, ...
            m.DNAinfo{i}{3}, m.DNAinfo{i}{4}, m.DNAinfo{i}{5}, 'Setup Reactions');
        m.DNAinfo{i}{6} = 'rxns_already_set_up';
    end
end % end of first run

if m.FirstRun
    txtl_enzyme_resource_degradation(modelObj);
    m.FirstRun = false;
end

set(modelObj, 'UserData', m)

%% SUBSEQUENT RUNS
if ~isempty(time_vector) && size(time_vector,1) > 1
    prevData = zeros(size(time_vector,1),size(modelObj.Species,1));
end

% Species-data pairs is needed
if iscell(data)
    SpName = findspecies(modelObj,data{:,1});
    for k=1:size(data,1)
        if size(data{k,2},1) == 1
            %first run initial amount provided
            modelObj.Species(SpName(k)).InitialAmount = normalizeSmallAmount(data{k,2});
        elseif size(data{k,2},1) > 1
            %setting up the initial amount to the latest simulation result
            modelObj.Species(k).InitialAmount = normalizeSmallAmount(data{k,2}(end));
            % reordering the data cell matrix to make it compatible with order in modelObj
            prevData(:,SpName(k)) = data{k,2}(:);
            
        else
            % if no data was given, we issue a warning
            warning('something went wrong on data(%d,2)',k);
        end
    end
    % we have as many colums of species in the model as in data, we set the
    % inital amount to that (here the order of data matters!)
    
elseif size(modelObj.Species,1) == size(data,2) && size(data,1) == 1
    for k=1:size(data,2)
        modelObj.Species(k).InitialAmount = normalizeSmallAmount(data(1,k));
    end
    % we have simulation data from a previous run
elseif size(modelObj.Species,1) == size(data,2) && size(data,1) > 1
    for k=1:size(data,2)
        modelObj.Species(k).InitialAmount = normalizeSmallAmount(data(end,k));
        prevData(:,k) = data(:,k);
    end
else
    % no data was provided, no action needed
end


% initial amounts set in modelObj.Species(k).InitialAmount.
% previousdata, if any, stored in prevData.
if ~isempty(configsetObj)
    
    simData = sbiosimulate(modelObj, configsetObj);
    
    if isempty(time_vector)
        [t_ode, prelimData, metaData] = getdata(simData, 'nummetadata');
        species_ind = cellfun(@(x) strcmp(x.Type,'species'), metaData);
        x_ode = prelimData(:,species_ind);
    else
%         %simData.Data also contains the time series for nonconstant
%         %parameters 
%         %(ie, it has 'AGTPdeg_F' data too, and this is a parameter, not a species). 
%         %So we need to explicitly pick out just the species.
%         newData = zeros(length(simData.Time), size(modelObj.species,1));
%         for i = 1:size(modelObj.species,1)
%         idx = strcmp(simData.DataNames, modelObj.species(i).Name);
%         newData(:,i) = simData.Data(:,idx);
%         end
        [new_Time, prelimData, metaData] = getdata(simData, 'nummetadata');
        species_ind = cellfun(@(x) strcmp(x.Type,'species'), metaData);
        new_Data = prelimData(:,species_ind);
        
        t_ode = [time_vector; new_Time+time_vector(end)];
        x_ode = [prevData;new_Data];
    end
    
    % TODO zoltuz 03/05/13 we need a new copyobject for simData -> merge
    % two simData object.
    varargout{1} = simData;
    
else
    
    % parameter estimation mode, no simulation
    x_ode = [];
    t_ode = [];
    simData = [];
end

switch nargout
    case 0
        varargout{1} = [];
    case 1
        varargout{1} = simData;
    case 2
        varargout{1} = t_ode;
        varargout{2} = x_ode;
    case 3
        varargout{1} = t_ode;
        varargout{2} = x_ode;
        varargout{3} = modelObj;
    case 4
        varargout{1} = t_ode;
        varargout{2} = x_ode;
        varargout{3} = modelObj;
        varargout{4} = simData;
    otherwise
        error('not supported operation mode');
        
end
end

% InitialAmount cannot be smaller than certain amount, therefore small
% amounts are converted to zero
function retValue = normalizeSmallAmount(inValue)

if (abs(inValue) < eps) || inValue<0 % treats values below eps as zero.
    retValue = 0;
else
    retValue = inValue;
end

end
