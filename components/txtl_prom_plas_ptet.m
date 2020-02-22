% txtl_prom_plas_ptet.m - promoter information for plasR and ptet combinatorial promoter
% Vipul Singhal Dec 2013
%
% This file contains a description of the plasR and ptet combinatorial promoter.
% Calling the function txtl_prom_plasR_ptet() will set up the reactions for
% transcription with the measured binding rates and transription rates.
% 
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

function varargout= txtl_prom_plas_ptet(mode, tube, dna, rna, varargin)


    % Create strings for reactants and products
    DNA = ['[' dna.Name ']'];		% DNA species name for reactions
    RNA = ['[' rna.Name ']'];		% RNA species name for reactions
    RNAP = 'RNAP';			% RNA polymerase name for reactions
    RNAPbound = ['RNAP:' dna.Name];
%     P2 = 'protein tetRdimer'; % we dont know a priori which versions of
%     the activator and repressor proteins will be present. 
%     P3 = 'OC12HSL:protein lasR';
%     DNAP3 = [ dna.Name ':' P3 ];
    % importing the corresponding parameters
    paramObj = txtl_component_config('plas_ptet');
    
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
    defaultBasePairs = {'plas_ptet','junk','thio';150,500,0};
    promoterData = txtl_setup_default_basepair_length(tube,promoterData,...
        defaultBasePairs);
    
    varargout{1} = promoterData;
    
    coreSpecies = {RNAP, RNAPbound};
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)),'Internal');

%     if mode.utr_attenuator_flag
%         txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' P3 ],prom_spec, rbs_spec, gene_spec,{P3} ); 
%     else
%         txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' P3 ],{P3}); 
%     end


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

    % leaky expression (put in a leaky expression mode)
    
    %find the activator species and the repressor species
    matchStr = regexp(listOfSpecies,'(^protein tetR.*dimer$)','tokens','once'); 
    listOftetRdimer = vertcat(matchStr{:});
    
    
    % June 8, 2019: add leaky expression. Maybe in a future release this
    % will be a mode. But for now we just have it as the default. The leaky
    % expression parameter was estimated to be 35 in log space. This is the
    % value being put into the txtlsim paper. 
    % Note that based on the logic of the combinatorial promoter, 
    % this parameter value cannot be the ptet promoter
    % unrepressed expression value. that is too high to be leaky. This is
    % the plas leaky espression value.
    %
    %   leaky expression
    parameters = {'TXTL_PLAS_RNAPbound_F',paramObj.RNAPbound_Forward;...
        'TXTL_PLAS_RNAPbound_R',paramObj.RNAPbound_Reverse};
    txtl_addreaction(tube,[DNA ' + ' RNAP ' <-> [' RNAPbound ']'],...
        'MassAction',parameters);
    
    p = regexp(listOfSpecies,'^OC12HSL:protein lasR(-lva)?$', 'match');
    listOfactivators = vertcat(p{:});
    
    % setup the activation and repression reactions. 
    % logic: repression supersedes activation. Even if we have an activated
    % promoter, if the repressor binds to it later, it will cause the
    % activator to fall off. 
    
    % first define the activator reactions, since the repression reactions
    % depend on these. again, will not be a problem in the multipass case. 
    for k = 1:size(listOfactivators,1)
        setup_promoter_reactions(mode, tube, dna, rna, RNAP,RNAPbound,...
            prom_spec, rbs_spec, gene_spec,listOfactivators{k}, paramObj)
    end
    
    % to pick the species to which we want tetR to be bound, we
    % search over the set of bound activator - DNA complexes, without tetR
    % bound
    listOfSpecies_justAdded = get(tube.species, 'name');
    activ_bound = ['^(RNAP:' dna.name '|' dna.name '):(OC12HSL:protein lasR(-lva)?)$'];
%     tetR dimer binds to the string after (at the end of the string) the activated lasR. and
%     the AGTP and CUTP bind before (ie, at the start of the string, / earlier in the string). 
%     I.e., dont want the following species to to get recognized
%     'RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR:protein tetRdimer'
%     'AGTP:CUTP:RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR:protein tetRdimer'
%     'AGTP:CUTP:RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
%     'AGTP:RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
%     but do want the following to bind: 
%     'DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
%     'RNAP:DNA plas_ptet--utr1--deGFP:OC12HSL:protein lasR'
    
    [tkns, species_list_2] = regexp(listOfSpecies_justAdded,activ_bound,'tokens', 'match');
    
    list_of_repressible_dna = vertcat(species_list_2{:});
    list_of_tokens = vertcat(tkns{:});
    
    % regexp logic is: start with dna with the plas promoter, or with RNAP
    % bound to the plas promoter. end with the activated protein bound. 
    % 
    matchStr = regexp(listOfSpecies_justAdded,'(^protein tetR.*dimer$)','tokens','once'); 
    listOftetRdimer = vertcat(matchStr{:});
    % repression of ptet by tetR dimer
    
    for k = 1:size(listOftetRdimer,1)
        
        % repressor binding to DNA
        txtl_addreaction(tube,...
            ['[' dna.name '] + ' listOftetRdimer{k} ' <-> [' dna.name ':' listOftetRdimer{k} ']'],...
            'MassAction',{'TXTL_PTET_sequestration_F',getDNASequestrationRates(paramObj,'F');...
            'TXTL_PTET_sequestration_R',getDNASequestrationRates(paramObj,'R')});
        for i = 1:size(list_of_repressible_dna)
            
            %repressor binding to activated dna
            txtl_addreaction(tube,...
                [list_of_repressible_dna{i} ' + ' listOftetRdimer{k} ' <-> ['...
                list_of_repressible_dna{i} ':' listOftetRdimer{k} ']'],...
                'MassAction',{'TXTL_PTET_sequestration_F',getDNASequestrationRates(paramObj,'F');...
                          'TXTL_PTET_sequestration_R',getDNASequestrationRates(paramObj,'R')});
                      
%            repressor knocking off activator
            % use the tokens captured previously to set up the reacion
            % products. (DNA, RNAP if present, and ligand bound activator)
            % first token is the RNAP optionally bound to the dna, and the
            % second is the activator complex. 
            
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
    
    
    
    % transcription only happens when there is no repression and there is
    % positive activation.
    
%     
%     %% plasR activaton
%     Robj3 = addreaction(tube, [dna.Name ' + ' P3 ' <-> [' dna.Name ':' P3 ']' ]);
%     Kobj3 = addkineticlaw(Robj3, 'MassAction');
%     addparameter(Kobj3, 'kf', paramObj.activation_F);
%     addparameter(Kobj3, 'kr', paramObj.activation_R);
%     set(Kobj3, 'ParameterVariableNames', {'kf', 'kr'});
%      
%     Robj6 = addreaction(tube, [dna.Name ':' P3 ' + ' RNAP ' <-> [' RNAPbound ':' P3 ']' ]);
%     Kobj6 = addkineticlaw(Robj6, 'MassAction');
%     addparameter(Kobj6, 'kf', paramObj.RNAPbound_Forward);
%     addparameter(Kobj6, 'kr', paramObj.RNAPbound_Reverse);
%     set(Kobj6, 'ParameterVariableNames', {'kf', 'kr'});
   
    %%
%     if mode.utr_attenuator_flag
%         txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP,...
%             [RNAPbound ':' P3],prom_spec, rbs_spec, gene_spec,{P3} );
%     else
%         txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' P3 ],{P3});  
%     end
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%   
else
    error('txtltoolbox:txtl_prom_plas_ptet:undefinedmode', ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end 
end

function setup_promoter_reactions(mode, tube, dna,rna, RNAP, RNAPbound,...
    prom_spec, rbs_spec, gene_spec, TF, paramObj)
% !TODO generalize this function to be used with all activator files. remove any
% plas specificity. VS aug 2017 -- Note to future self: an activator class
% would work really well here. Basically we must have activator, repressor,
% etc classes. 

txtl_addreaction(tube, ...
    [dna.Name ' + ' TF ' <-> [' dna.Name ':' TF ']' ],...
    'MassAction',{'TXTL_PLAS_TFBIND_F',paramObj.activation_F;...
    'TXTL_PLAS_TFBIND_R',paramObj.activation_R});

txtl_addreaction(tube, ...
    [dna.Name ':' TF ' + ' RNAP ' <-> [' RNAPbound ':' TF ']' ],...
    'MassAction',{'TXTL_PLAS_TFRNAPbound_F',paramObj.RNAPbound_Forward_actv;...
    'TXTL_PLAS_TFRNAPbound_R',paramObj.RNAPbound_Reverse_actv});

txtl_addreaction(tube, ...
    [RNAPbound ' + ' TF ' <-> [' RNAPbound ':' TF ']' ],...
    'MassAction',{'TXTL_PLAS_TFBIND_F',paramObj.activation_F;...
    'TXTL_PLAS_TFBIND_R',paramObj.activation_R}); % this line was added Jun 8, 2019
% 
% if mode.utr_attenuator_flag
%     mode.add_dna_driver = 'Setup Species';
%     txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' TF],prom_spec, rbs_spec, gene_spec,{TF} );
%     mode.add_dna_driver = 'Setup Reactions';
%     txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' TF],prom_spec, rbs_spec, gene_spec,{TF} );
% else
    mode.add_dna_driver = 'Setup Species';
    
    txtl_transcription(mode, tube, dna, rna, RNAP,RNAPbound);  %added on Jun 30, 2019
    txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' TF ],{TF});
    mode.add_dna_driver = 'Setup Reactions';
    
    txtl_transcription(mode, tube, dna, rna, RNAP,RNAPbound); %added on Jun 30, 2019
    txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' TF ],{TF});
% end

end


% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
