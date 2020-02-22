% Written by Vipul Singhal Nov 2012
% 
% Copyright (c) 2012 by California Institute of Technology
% All rights reserved.
%
% Documentation:
% there are 4 ways to call this function. 
% [t_ode_output, x_ode_output, simData_output, t_sen_output, x_sen_output, senOutputs, senInputs] = txtl_runsim(modelObj, configsetObj, t_ode_input, x_ode_input, t_sen_input, x_sen_input)
% [t_ode_output, x_ode_output, simData_output, t_sen_output, x_sen_output, senOutputs, senInputs] = txtl_runsim(modelObj, configsetObj, t_ode_input, x_ode_input, t_sen_input, x_sen_input, simData_input)
% where the first 3 outputs are mandatory

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

function [t_ode, x_ode, simData, varargout] = txtl_runsim_sensitivity(varargin)
        

    if nargin ~= 6 && nargin~= 7
        error('txtl_runsim_sensitivity should be called either with 6 or 7 arguments!');
    else
        modelObj = varargin{1};
        configsetObj = varargin{2};
        time_vector = varargin{3};
        data = varargin{4};
        t_sen_input = varargin{5};
        x_sen_input = varargin{6};
        if nargin == 7
            simData = varargin{7};
        end
        if configsetObj.SolverOptions.SensitivityAnalysis ~= 1
            % in case runsim is called in sensitivityAnalysis mode but without 
            % the sensitivity analysis configset. 
            configsetObj.SolverOptions.SensitivityAnalysis = 1;
        end        
    end

if ~isempty(time_vector) && size(time_vector,1) > 1
    prevData = zeros(size(time_vector,1),size(modelObj.Species,1));
end
% Species-data pairs is needed
if iscell(data) 
    SpName = findspecies(modelObj,data{:,1}); 
    for k=1:size(data,1)
       if size(data{k,2},1) == 1
           %first run initial amount provided
           modelObj.Species(k).InitialAmount = normalizeSmallAmount(data{k,2});
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


% run the simbiology simulation

simData = sbiosimulate(modelObj, configsetObj);

% get the indices for the sensitivity data and the species evolution data
h = get(simData, 'DataInfo');
speciesCount = 1;
sensitivityCount = 1;
for i = 1:length(h)
    if strcmp(h{i}.Type,'species')
        speciesIndex(speciesCount) = i;
        speciesCount = speciesCount+1;
         
    elseif strcmp(h{i}.Type,'sensitivity')
        sensitivityIndex(sensitivityCount) = i;
        sensitivityCount = sensitivityCount+1;
    end
end
    

if isempty(time_vector) % if there is no previous data to append
    x_ode = simData.Data;
    t_ode = simData.Time;
    [t_sen,x_sen, sensOutputs, sensInputs] = getsensmatrix(simData);

else % if data from previous simulations needs to be appended
     % put in a check for whether the data from prev runs and the current runs have the same indices. 
     
        t_ode = [time_vector; simData.Time+time_vector(end)];
        x_ode = [prevData;simData.Data(:,speciesIndex)]; 
        [t_sen_new,x_sen_new, sensOutputs, sensInputs] = getsensmatrix(simData);
        %display('The size of the sensitivity data')
        size(t_sen_new)
        
        size(t_sen_input)
       
        t_sen = [t_sen_input; t_sen_input(end)+t_sen_new];
        x_sen = [x_sen_input; x_sen_new];
end

    varargout{1} = t_sen;
    varargout{2} = x_sen;
    varargout{3} = sensOutputs;
    varargout{4} = sensInputs;

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