% txtl_prom_pBAD_ptet.m - promoter information for pBAD and ptet combinatorial promoter
% Zoltan Tuza, Oct 2012
% Vipul Singhal Jun 2014, Aug 2017
% 
%
% This file contains a description of the pBAD and ptet combinatorial promoter.
% Calling the function txtl_prom_pBAD_ptet() will set up the reactions for
% transcription with the measured binding rates and transription rates.
% This file is based on the file txtl_prom_plas_ptet.m
% 

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

function varargout= txtl_prom_pBAD_ptet(mode, tube, dna, rna, varargin)


    % Create strings for reactants and products
    DNA = ['[' dna.Name ']'];		% DNA species name for reactions
    RNA = ['[' rna.Name ']'];		% RNA species name for reactions
    RNAP = 'RNAP';			% RNA polymerase name for reactions
    RNAPbound = ['RNAP:' dna.Name];
%     P1 = 'protein sigma70';
    
%     P2 = 'protein tetRdimer';
%     P3 = 'protein AraC';
%     AraCbound = ['arabinose:' P3];
    
    % importing the corresponding parameters
    paramObj = txtl_component_config('pBAD_ptet');
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode.add_dna_driver, 'Setup Species')
    
    promoterData = varargin{1};
   if nargin==8
    prom_spec = varargin{2};
    rbs_spec = varargin{3};
    gene_spec = varargin{4};
    elseif nargin~=5
        error('the number of argument should be 5 or 8, not %d',nargin);
    end
    defaultBasePairs = {'pBAD_ptet','junk','thio';150,500,0};
    promoterData = txtl_setup_default_basepair_length(tube,promoterData,...
        defaultBasePairs);
    
    varargout{1} = promoterData;
    
    coreSpecies = {RNAP,RNAPbound};

    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)),'Internal');
    % empty cellarray for amount => zero amount
%     if mode.utr_attenuator_flag
%         txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, RNAPbound, prom_spec, rbs_spec, gene_spec );
%         txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' P2 ], prom_spec, rbs_spec, {P2} );
%         txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' AraCbound ], prom_spec, rbs_spec, gene_spec,{AraCbound} );
%         txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' P2 ':' AraCbound ], prom_spec, rbs_spec, gene_spec,{P2, AraCbound} );
%     else
%         txtl_transcription(mode, tube, dna, rna, RNAP,RNAPbound); %leaky slow rate
%         txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' P2 ],{P2}); %lowest rate
%         txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' AraCbound ],{AraCbound}); %highest rate
%         txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' P2 ':' AraCbound ],{P2, AraCbound}); %slightly higher than 1.
%     end

    %(check agains shaobin results. the parameters here should be tuned to
    %get the shaobin curves. translation/degradation etc should be standard.)

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif strcmp(mode.add_dna_driver,'Setup Reactions')
    listOfSpecies = varargin{1};
    if nargin==8
    prom_spec = varargin{2};
    rbs_spec = varargin{3};
    gene_spec = varargin{4};
    elseif nargin~=5
        error('the number of argument should be 5 or 8, not %d',nargin);
    end
    
    matchStr = regexp(listOfSpecies,'(^protein tetR.*dimer$)','tokens','once'); 
    listOftetRdimer = vertcat(matchStr{:});
    
    p = regexp(listOfSpecies,'^arabinose:protein AraC(-lva)?$', 'match');
    listOfactivators = vertcat(p{:});    
    
    for k = 1:size(listOfactivators,1)
        setup_promoter_reactions(mode, tube, dna, rna, RNAP,RNAPbound,...
            prom_spec, rbs_spec, gene_spec,listOfactivators{k}, paramObj)
    end    
    
    listOfSpecies_justAdded = get(tube.species, 'name');
    activ_bound = ['^(RNAP:' dna.name '|' dna.name '):(arabinose:protein AraC(-lva)?)$'];
    
    [tkns, species_list_2] = regexp(listOfSpecies_justAdded,activ_bound,'tokens', 'match');
    
    list_of_repressible_dna = vertcat(species_list_2{:});
    list_of_tokens = vertcat(tkns{:});    
    
    matchStr = regexp(listOfSpecies_justAdded,'(^protein tetR.*dimer$)','tokens','once'); 
    listOftetRdimer = vertcat(matchStr{:});
    % repression of ptet by tetR dimer
    
    for k = 1:size(listOftetRdimer,1)
        
        % repressor binding to DNA
        txtl_addreaction(tube,...
            ['[' dna.name '] + ' listOftetRdimer{k} ' <-> [' dna.name ':' listOftetRdimer{k} ']'],...
            'MassAction',{'ptet_sequestration_F',getDNASequestrationRates(paramObj,'F');...
            'ptet_sequestration_R',getDNASequestrationRates(paramObj,'R')});
        for i = 1:size(list_of_repressible_dna)
            
            %repressor binding to activated dna
            txtl_addreaction(tube,...
                [list_of_repressible_dna{i} ' + ' listOftetRdimer{k} ' <-> ['...
                list_of_repressible_dna{i} ':' listOftetRdimer{k} ']'],...
                'MassAction',{'ptet_sequestration_F',getDNASequestrationRates(paramObj,'F');...
                          'ptet_sequestration_R',getDNASequestrationRates(paramObj,'R')});
                      
%            repressor knocking off activator
            % check if RNAP is present in the first token
            token1 = list_of_tokens{i}{1};
            token2 = list_of_tokens{i}{2};
            if isempty(regexp(token1, 'RNAP', 'once'))
                txtl_addreaction(tube,...
                    ['[' list_of_repressible_dna{i} ':' listOftetRdimer{k} '] <-> ['...
                    dna.name ':' listOftetRdimer{k} '] + [' token2 ']'],...
                    'MassAction',{'TXTL_COMBINATORIAL_ACTIVATOR_KNOCKOFF_F',paramObj.Combinatorial_activ_knockoff_F;...
                                'TXTL_COMBINATORIAL_ACTIVATOR_KNOCKOFF_R',paramObj.Combinatorial_activ_knockoff_R});   
            else
                txtl_addreaction(tube,...
                    ['[' list_of_repressible_dna{i} ':' listOftetRdimer{k} '] <-> ['...
                    dna.name ':' listOftetRdimer{k} '] + RNAP + [' token2 ']'],...
                    'MassAction',{'TXTL_COMBINATORIAL_ACTIVATOR_KNOCKOFF_F',paramObj.Combinatorial_activ_knockoff_F;...
                                'TXTL_COMBINATORIAL_ACTIVATOR_KNOCKOFF_R',paramObj.Combinatorial_activ_knockoff_R});                  
            end
        end
    end
    
    
    
    
    
    
    
    
    
    
    
    %%
%     % parameters for leaky transcription
%     parameters = {'TXTL_PBADPTET_RNAPbound_F',paramObj.RNAPbound_Forward;...
%                   'TXTL_PBADPTET_RNAPbound_R',paramObj.RNAPbound_Reverse};
%     % Set up binding reaction
%     txtl_addreaction(tube,[DNA ' + ' RNAP ' <-> [' RNAPbound ']'],...
%         'MassAction',parameters);
    %

%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%   
else
    error('txtltoolbox:txtl_prom_pBAD_ptet:undefinedmode', ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end
end


function setup_promoter_reactions(mode, tube, dna,rna, RNAP, RNAPbound,...
    prom_spec, rbs_spec, gene_spec, TF, paramObj)
% !TODO generalize this function to be used with all activator files. remove any
% plas specificity. VS aug 2017

txtl_addreaction(tube, ...
    [dna.Name ' + ' TF ' <-> [' dna.Name ':' TF ']' ],...
    'MassAction',{'TXTL_PBAD_TFBIND_F',paramObj.activation_F;...
    'TXTL_PBAD_TFBIND_R',paramObj.activation_R});

txtl_addreaction(tube, ...
    [dna.Name ':' TF ' + ' RNAP ' <-> [' RNAPbound ':' TF ']' ],...
    'MassAction',{'TXTL_PBAD_TFRNAPbound_F',paramObj.RNAPbound_Forward_actv;...
    'TXTL_PBAD_TFRNAPbound_R',paramObj.RNAPbound_Reverse_actv});

% 
% if mode.utr_attenuator_flag
%     mode.add_dna_driver = 'Setup Species';
%     txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' TF],prom_spec, rbs_spec, gene_spec,{TF} );
%     mode.add_dna_driver = 'Setup Reactions';
%     txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' TF],prom_spec, rbs_spec, gene_spec,{TF} );
% else
    mode.add_dna_driver = 'Setup Species';
    txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' TF ],{TF});
    mode.add_dna_driver = 'Setup Reactions';
    txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' TF ],{TF});
% end

end

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
