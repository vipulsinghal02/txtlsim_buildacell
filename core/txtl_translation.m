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


function txtl_translation(mode, tube, dna, rna, protein, Ribobound)
% basically, this function does nothing if the RNA does not have a ribosome
% binding site. so for instance if the RNA is: [RNA att] or [RNA anti].


%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode.add_dna_driver, 'Setup Species')
    % Set up the species for translation
    Ribobound_term = ['term_' Ribobound.Name ];
    coreSpecies = {'AA',['AA:AGTP:' Ribobound.Name],Ribobound_term, 'Ribo'};
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');
    
    
    %%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode.add_dna_driver, 'Setup Reactions')
    
    %% Resource binding
    AAparameters = {'TL_AA_F',tube.UserData.ReactionConfig.TL_AA_Forward;
        'TL_AA_R',tube.UserData.ReactionConfig.TL_AA_Reverse};
    AGTPparameters = {'TL_AGTP_F',tube.UserData.ReactionConfig.TL_AGTP_Forward;
        'TL_AGTP_R',tube.UserData.ReactionConfig.TL_AGTP_Reverse};
    
    % resource binding
    txtl_addreaction(tube, ...
        ['[' Ribobound.Name '] + AA <-> [AA:' Ribobound.Name ']'],...
        'MassAction',AAparameters);
    txtl_addreaction(tube, ...
        ['[AA:' Ribobound.Name ']  + AGTP <-> [AA:AGTP:' Ribobound.Name ']'],...
        'MassAction',AGTPparameters);
    
    
    %% create the rules for the global parameters for the TL reaction and the consumption reaction
    
    % add global elongation parameter
    tlglob = sbioselect(tube, 'Name', 'TL_elong_glob', 'Type', 'Parameter');
    if isempty(tlglob)
        addparameter(tube, 'TL_elong_glob',tube.UserData.ReactionConfig.Translation_Rate);
    end
    
    % grab the name strings for the protein
    temp = regexp(protein.name, 'protein (\w*)', 'tokens');
    protspec = [temp{1}{1}];
    
    % use this sting to name the translation and consumption reactions
    tlparamname = ['TL_translation_' protspec];
    resourceconsname = ['TL_REScons_' protspec];
    
    % get the protein length, which decides what the actual protein production
    % rate is
    proteinlength = round(protein.UserData); % this is in amino acids not nucleotides.
    
    % add the transcription parameter in the model scope, with the length
    % adjusted value. The value specified here actually does not matter
    % because we will use a rule to set it.
    addparameter(tube, tlparamname,0);
    
    % add the aa consumption parameter in the global scope. the value it is
    % initialized to does not matter, since a rule will set it.
    addparameter(tube, resourceconsname, 0);
    
    % now we actually add the rule that sets the translation rate and the
    % aa consumption rate.
    ruleStr = [tlparamname ' =  TL_elong_glob/' num2str(proteinlength)];
    if isempty(sbioselect(tube,'Type','Rule', 'Rule', ruleStr))
        addrule(tube, ruleStr, 'initialAssignment');
    end
    
    % do the same for the consumption reactions
    ruleStr = [resourceconsname ...
        ' =  TL_elong_glob/' num2str(proteinlength) '*(' num2str(proteinlength) '-1)'];
    if isempty(sbioselect(tube,'Type','Rule', 'Rule', ruleStr))
        addrule(tube, ruleStr, 'initialAssignment');
    end
    
    
    Ribobound_term = ['term_' Ribobound.Name ];
    % add the reaction
    if isfield(tube.UserData, 'energymode') && strcmp(tube.UserData.energymode, 'regeneration')
        
        reactionObj = addreaction(tube, ...
            ['[AA:AGTP:' Ribobound.Name '] -> ' Ribobound_term ' + ' protein.Name ' +  AGMP' ]);
        addkineticlaw (reactionObj, 'MassAction');
        reactionObj.KineticLaw.ParameterVariableNames = tlparamname;
        
        % add the consumption reactions.
        reactionObj = addreaction(tube, ...
            ['[AA:AGTP:' Ribobound.Name '] -> ' Ribobound_term ' +  AGMP']);
        addkineticlaw (reactionObj, 'MassAction');
        reactionObj.KineticLaw.ParameterVariableNames = resourceconsname;
        
        
    else
        reactionObj = addreaction(tube, ...
            ['[AA:AGTP:' Ribobound.Name '] -> ' Ribobound_term ' + ' protein.Name ]);
        addkineticlaw (reactionObj, 'MassAction');
        reactionObj.KineticLaw.ParameterVariableNames = tlparamname;
        
        % add the consumption reactions.
        reactionObj = addreaction(tube, ...
            ['[AA:AGTP:' Ribobound.Name '] -> ' Ribobound_term]);
        addkineticlaw (reactionObj, 'MassAction');
        reactionObj.KineticLaw.ParameterVariableNames = resourceconsname;
    end
    
    
    %%%%%
    
    
    % translation termination reaction
    txtl_addreaction(tube,['[' Ribobound_term '] -> ' rna.Name ' +  Ribo'],...
        'MassAction',{'TXTL_RIBOBOUND_TERMINATION_RATE', tube.UserData.ReactionConfig.Ribobound_termination_rate});
    %%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error('txtltoolbox:txtl_translation:undefinedmode', ...
        'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end


end