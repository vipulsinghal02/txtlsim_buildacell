%TXTL_COMBINE	Generate a model by composition of two or more submodels
%
% newtube = txtl_combine(tubelist, vollist) - combine the contents of
% a list of tubes into a new tube.  The volume taken from each tube is
% given in vollist.  All species, reactions, events and other parameters
% are combined in the new tube.

% Written by Richard Murray, Sep 2012
%
% Copyright (c) 2012 by California Institute of Technology
% All rights reserved.
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

function Mobj = txtl_combine(tubelist, varargin)

if nargin == 2
    warning('You are changing the Extract-Buffer-DNA ratios, is that what you really want?'); 
    vollist = varargin{1};
else
    vollist = [10/3 4.41666 2.25];
end
    

% Create a model for the mixture
Mobj = txtl_newtube('mix_of_');

for m = 1:length(tubelist)
  tube = tubelist(m);
  Mobj.Name = [Mobj.Name tube.Name '_'];
  % Start with species, keeping track of total amount
  species_list = get(tube, 'Species');
  for i = 1:length(species_list)
    % Check first to see if this species already exists
    name = get(species_list(i), 'Name');
    species = getspecies(Mobj, name);
    if (isempty(species))
      % Create a new species in this model
      species = copyobj(species_list(i), Mobj);
      
      % Adjust amount by volume (re-normalized at bottom)
      species.InitialAmount = species.InitialAmount * vollist(m);
    else
      % Increase the concentration, scaled by volume (normalized at bottom)
      species.InitialAmount = species.InitialAmount + ...
	species_list(i).InitialAmount * vollist(m);
      if ~isempty(species.UserData)
          warning('%s.UserData was not empty!',species.Name)
      end
      species.UserData = species_list(i).UserData;
    end
  end

  % Copy the rest of the contents of the model
  copyobjlist(tube.Reactions, Mobj);
  copyobjlist(tube.Parameters, Mobj);
  copyobjlist(tube.Events, Mobj);
  copyobjlist(tube.Rules, Mobj);
  
  mobjUser = get(Mobj, 'UserData');
  tubeUser = get(tube, 'UserData');  
  if ~isempty(tubeUser.ReactionConfig)%, 'class'
      if  ~isempty(mobjUser.ReactionConfig)%, 'class'
          warning('Warning:ReactionConfigAlreadyPresent','Reaction Config field already populated, will be overwritten.')
          mobjUser.ReactionConfig = tubeUser.ReactionConfig;
      else
          mobjUser.ReactionConfig = tubeUser.ReactionConfig;
      end
  end
  if ~isempty(tubeUser.DNAinfo)
      tubeDNA = tubeUser.DNAinfo;
      mobjDNA = mobjUser.DNAinfo;
      if ~isempty(mobjDNA)
        outDNA = {mobjDNA(:); tubeDNA(:)};
      else
          outDNA = tubeDNA;
      end
      mobjUser.DNAinfo = outDNA;
  end
  mobjUser.FirstRun = true;
  set(Mobj,'UserData', mobjUser);
%   if ~isempty(tubeUser)
%       if isempty(mobjUser)
%           set(Mobj, 'UserData', tubeUser);
%       else
%           if isa(mobjUser, 'txtl_reaction_config')
%               outCell = {mobjUser(:); tubeUser(:)};
%               set(Mobj, 'UserData', outCell);
% %           elseif iscell(mobjUser)
% %               outCell = mobjUser{1};
% %               for k = 1:length(mobjUser)
% %                   outCell = {outcell mobjUser{:}; tubeUser(:)};
% %                   set(Mobj, 'UserData', outCell);
% %               end
%           end
%       end
%   end
end

Mobj.Name = Mobj.Name(1:end-1); % deleting the last '_' character
% Go through and normalize species concentration
totalvol = sum(vollist(1:length(tubelist)), 'double');
for i = 1:length(Mobj.Species)
  Mobj.Species(i).InitialAmount = Mobj.Species(i).InitialAmount / totalvol;
end

Mobj.UserData.energymode = 'regeneration';  
% See:  
% https://www.mathworks.com/help/simbio/ug/selecting-absolute-  
% tolerance-and-relative-tolerance-for-simulation.html  
cs1 = getconfigset(Mobj); 
set(cs1.RuntimeOptions, 'StatesToLog', 'all');  
% set(cs1.SolverOptions, 'AbsoluteToleranceScaling', 1);  
% set(cs1.SolverOptions, 'AbsoluteTolerance', 1.0e-6);  
% set(cs1.SolverOptions, 'AbsoluteToleranceStepSize', 6*3600*1.0e-6*0.1); 
% set(cs1.SolverOptions, 'RelativeTolerance', 1.0e-6);  
% try: AbsoluteToleranceStepSize = AbsoluteTolerance * StopTime * 0.1 


return

% Copy a list of objects to a model
function copyobjlist(list, model)
for i = 1:length(list)
  copyobj(list(i), model);
end

% Search for a species in an existing model
function species = getspecies(model, name)
species_list = get(model, 'Species');
for i = 1:length(species_list)
  if strcmp(get(species_list(i), 'Name'), name) == 1
    species = species_list(i);
    return 
  end
end

species = [];
return
