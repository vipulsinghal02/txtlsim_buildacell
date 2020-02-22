% txtl_transcription.m - sigma factor independent implementation for gene
% transcription in the TXTL system
% RMM, 9 Sep 2012
%
% It can be called by promoter files that need to
% set up the approriate transcription reactions.

% Written by Richard Murray, 9 Sep 2012
% Edited by Vipul Singhal, 2012 - 2017
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

function txtl_transcription(mode, varargin)
tube = varargin{1};
dna = varargin{2};
rna = varargin{3};
RNAP = varargin{4}; % RNA polymerase name for reactions
RNAPbound = varargin{5};
% if called with extra species, then RNAPbound is more than just RNAP70:DNA. there are other sopecies also attached.

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode.add_dna_driver, 'Setup Species')
    RNAPbound_term = ['term_' RNAPbound];
    if nargin < 6
        error('the number of argument should be at least 6, not %d',nargin);
    elseif nargin > 6
        extraSpecies = varargin{6};
        coreSpecies = {'CUTP', 'AGTP',RNAPbound,['CUTP:AGTP:' RNAPbound],RNAPbound_term,RNAP,extraSpecies{:}};
    else
        coreSpecies = {'CUTP', 'AGTP',RNAPbound,['CUTP:AGTP:' RNAPbound],RNAPbound_term,RNAP};
    end
    
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');
    
    %%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode.add_dna_driver,'Setup Reactions')
    
    % calculate the transcription rate from information in the config file
    % and the length of the gene to be transcribed
    
    % add global elongation parameter,
    txglob = sbioselect(tube, 'Name', 'TX_elong_glob', 'Type', 'Parameter');
    if isempty(txglob)
        addparameter(tube, 'TX_elong_glob',tube.Userdata.ReactionConfig.Transcription_Rate);
    end
    
    
    % then add the dependent parameters for both the tx and the consumption
    % reactions, also at the global scope, and tie the three parameters via
    % a couple of initial assignment rules.
    
    % start with grabbing the name strings for the RNA
    temp = regexp(rna.name, 'RNA (\w*)--(\w*)', 'tokens');
    rnaspec = [temp{1}{1} '_' temp{1}{2}];
    
    % use this sting to name the transcription and consumption reactions
    txparamname = ['TX_transcription_' rnaspec];
    ntpconsname = ['TX_NTPcons_' rnaspec];
    
    % get the RNA length, which decides what the actual mRNA production
    % rate is
    RNAlength = rna.UserData;
    % add the transcription parameter in the model scope, with the length
    % adjusted value.
    if isempty(sbioselect(tube, 'Type', 'Parameter', 'Name', txparamname))
        addparameter(tube, txparamname,tube.Userdata.ReactionConfig.Transcription_Rate/RNAlength);
    end
    
    % compute the consumption reaction rate as follows
    % ntpcnt = rna.length/4.
    ntpcnt = round(RNAlength/4); %
    % add the ntp consumption parameter in the global scope.
    if isempty(sbioselect(tube, 'Type', 'Parameter', 'Name', ntpconsname))
        addparameter(tube, ntpconsname,tube.Userdata.ReactionConfig.Transcription_Rate/RNAlength*(ntpcnt-1));
    end
    
    
    % how to think about this: 1 nM of AGTP is 0.5nM of ATP, 0.5nM of GTP.
    % Since the reaction rnap:dna:agtp:cutp -> rnap:dna_term + mrna uses
    % 2nM of NTPs to make 1nM of mRNA, if the mrna is 1000ntp long, then we
    % still need to consume 998nM of ntp. The rate at which the above
    % reaction took place is kt/1000.
    %
    % so now we consider the ntp consumption reaction.
    % rnap:dna:agtp:cutp -> rnap:dna
    % each time this reaction fires, it uses 2nM of ntp. so it needs to
    % fire 998/2 = 499 times for each time the earlier reaction fires. Now
    % 499 is 1000/2 - 1, thus we define ntpcnt = rnalength/2, and the
    % reaction of the consumption reaction is ktx/rnalength*(ntpcnt-1)
    
    % now we actually add the rule that sets the transcription rate and the
    % ntp consumption rate.
    ruleStr = [txparamname ' = TX_elong_glob/' num2str(RNAlength)];
    if isempty(sbioselect(tube,'Type','Rule', 'Rule', ruleStr))
        addrule(tube, ruleStr, 'initialAssignment');
    end
    
    ruleStr = [ntpconsname ...
        ' = TX_elong_glob/' num2str(RNAlength) '*(' num2str(ntpcnt) '-1)'];
    if isempty(sbioselect(tube,'Type','Rule', 'Rule', ruleStr))
        addrule(tube, ruleStr, 'initialAssignment');
    end
    
    % write down the string for the transcription equation depending on the
    % mode of transcription. ie, if we have the energy regeneration mode,
    % then we model AGMP, otherwise we do the usual... except, with and
    % without energy mode ends up being exactly the same for
    % teranscription. ...
    if isfield(tube.UserData, 'energymode') && strcmp(tube.UserData.energymode, 'regeneration')
        
            RNAPbound_term = ['term_' RNAPbound];
            transcriptionEq = ...
                ['[CUTP:AGTP:' RNAPbound '] -> ' RNAPbound_term ' + '...
                rna.Name];
            
    else
        RNAPbound_term = ['term_' RNAPbound];
        transcriptionEq = ...
            ['[CUTP:AGTP:' RNAPbound '] -> ' RNAPbound_term ' + ' rna.Name];
        
    end
    
    % add the actual transcription reaction. Note that we use addreaction (
    % a simbiology function) as opposed to the txtl toolbox wrapper
    % txtl_addreaction, because we want to specify a model scoped parameter
    % as a parmeter, and not create a reaction scoped parameter.
    reactionObj = addreaction(tube,transcriptionEq);
    addkineticlaw (reactionObj, 'MassAction');
    reactionObj.KineticLaw.ParameterVariableNames = txparamname;
    
    
    % add the consumption reactions
    reactionObj = addreaction(tube,['[CUTP:AGTP:' RNAPbound '] -> ' RNAPbound]);
    addkineticlaw (reactionObj, 'MassAction');
    reactionObj.KineticLaw.ParameterVariableNames = ntpconsname;
    
    
    
    % add the termination reactions
    if nargin < 6
        error('the number of argument should be at least 6, not %d',nargin);
    elseif nargin > 6
        extraSpecies = varargin{6};
        % processing the extraSpecies, like activators, inducers etc.
        extraStr = extraSpecies{1};
        for k=2:size(extraSpecies,2)
            extraStr = [extraStr ' + ' extraSpecies{k}];
        end
        
        % the termination reaction parameter is reaction scoped, and can be
        % globalized with the globalize_params function.
        txtl_addreaction(tube,['[' RNAPbound_term '] -> '...
            RNAP  ' + ' dna.Name ' + ' extraStr],...
            'MassAction',...
            {'TXTL_RNAPBOUND_TERMINATION_RATE', ...
            tube.UserData.ReactionConfig.RNAPbound_termination_rate});
        
    else
        %termination reaction
        txtl_addreaction(tube,['[' RNAPbound_term '] -> ' RNAP  ' + ' dna.Name],...
            'MassAction',...
            {'TXTL_RNAPBOUND_TERMINATION_RATE', ...
            tube.UserData.ReactionConfig.RNAPbound_termination_rate});
    end
    
    % define the nucleotide binding parameters
    NTPparameters_step1 = {'TXTL_NTP_RNAP_1_F', tube.UserData.ReactionConfig.NTP_Forward_1;
        'TXTL_NTP_RNAP_1_R', tube.UserData.ReactionConfig.NTP_Reverse_1};
    NTPparameters_step2 = {'TXTL_NTP_RNAP_2_F', tube.UserData.ReactionConfig.NTP_Forward_2;
        'TXTL_NTP_RNAP_2_R', tube.UserData.ReactionConfig.NTP_Reverse_2};
    
    % add the nucleotide binding reaction
    txtl_addreaction(tube,['[' RNAPbound '] + AGTP <-> [AGTP:' RNAPbound ']'],...
        'MassAction',NTPparameters_step1);
    txtl_addreaction(tube,['[' RNAPbound '] + CUTP <-> [CUTP:' RNAPbound ']'],...
        'MassAction',NTPparameters_step1);
    txtl_addreaction(tube,['[AGTP:' RNAPbound '] + CUTP <-> [CUTP:AGTP:' RNAPbound ']'],...
        'MassAction',NTPparameters_step2);
    txtl_addreaction(tube,['[CUTP:' RNAPbound '] + AGTP <-> [CUTP:AGTP:' RNAPbound ']'],...
        'MassAction',NTPparameters_step2);
    
    
    
    
    %%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error('txtltoolbox:txtl_transcription:undefinedmode', ...
        'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end




% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
% Parameters describing the enzymatic process
