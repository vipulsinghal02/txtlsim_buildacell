% txtl_prom_p28_ptet.m - promoter information for p28 and ptet combinatorial promoter
% Zoltan Tuza, Oct 2012
% Vipul Singhal Jun 2014
%
% This file contains a description of the p28 and ptet combinatorial promoter.
% Calling the function txtl_prom_p28_ptet() will set up the reactions for
% transcription with the measured binding rates and transription rates.
% 
% 

% Written by Zoltan Tuza, Oct 2012
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

function varargout= txtl_prom_p28_ptet(mode, tube, dna, rna, varargin)


    % Create strings for reactants and products
    DNA = ['[' dna.Name ']'];		% DNA species name for reactions
    RNA = ['[' rna.Name ']'];		% RNA species name for reactions
    RNAP = 'RNAP28';			% RNA polymerase name for reactions
    RNAPbound = ['RNAP28:' dna.Name];
    P1 = 'protein sigma28';
    P2 = 'protein tetRdimer';
    
    % importing the corresponding parameters
    paramObj = txtl_component_config('p28_tetR');
    
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
    defaultBasePairs = {'p28_ptet','junk','thio';150,500,0};
    promoterData = txtl_setup_default_basepair_length(tube,promoterData,...
        defaultBasePairs);
    
    varargout{1} = promoterData;
    
    coreSpecies = {RNAP,RNAPbound,P1,P2};
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');
    if mode.utr_attenuator_flag
        txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP,RNAPbound, prom_spec, rbs_spec, gene_spec ); %leaky slow rate
        txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' P2 ],prom_spec, rbs_spec, gene_spec,{P2} );  %lowest rate
    else
        txtl_transcription(mode, tube, dna, rna, RNAP,RNAPbound); %leaky slow rate
        txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' P2 ],{P2});  %lowest rate
    end

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif strcmp(mode.add_dna_driver,'Setup Reactions')
    
    % Parameters that describe this promoter
    parameters = {'TXTL_PTET_RNAPbound_F',paramObj.RNAPbound_Forward;...
                  'TXTL_PTET_RNAPbound_R',paramObj.RNAPbound_Reverse};
    % Set up binding reaction
    txtl_addreaction(tube,[DNA ' + ' RNAP ' <-> [' RNAPbound ']'],...
        'MassAction',parameters);

    % 
    % Set up bindig reaction for both sigma28 and tetR
    Robj3 = addreaction(tube, [dna.Name ' + ' P2 ' <-> ' dna.Name ':' P2 ]);
    Kobj3 = addkineticlaw(Robj3, 'MassAction');
    Pobj3f = addparameter(Kobj3, 'kf', 2.86e-3);
    Pobj3r = addparameter(Kobj3, 'kr', 0.11e-4);
    set(Kobj3, 'ParameterVariableNames', {'kf', 'kr'});
    % 
    Robj4 = addreaction(tube, [RNAPbound  ' + ' P2 ' <-> ' RNAPbound ':' P2 ]);
    Kobj4 = addkineticlaw(Robj4, 'MassAction');
    Pobj4f = addparameter(Kobj4, 'kf', 2.86e-3);
    Pobj4r = addparameter(Kobj4, 'kr', 0.11e-4);
    set(Kobj4, 'ParameterVariableNames', {'kf', 'kr'});
    % 
    % Set up binding reaction for tetR
    Robj2 = addreaction(tube, [dna.Name ':' P2 ' + ' RNAP ' <-> [' RNAPbound ':' P2 ']' ]);
    Kobj2 = addkineticlaw(Robj2, 'MassAction');
    Pobj2f = addparameter(Kobj2, 'kf', paramObj.RNAPbound_Forward);
    Pobj2r = addparameter(Kobj2, 'kr', paramObj.RNAPbound_Reverse*1000);
    set(Kobj2, 'ParameterVariableNames', {'kf', 'kr'});

    if mode.utr_attenuator_flag
        txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP,RNAPbound, prom_spec, rbs_spec, gene_spec ); %leaky slow rate
        txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' P2 ],prom_spec, rbs_spec, gene_spec,{P2} );  %lowest rate
    else
        txtl_transcription(mode, tube, dna, rna, RNAP,RNAPbound); %leaky slow rate
        txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' P2 ],{P2});  %lowest rate
    end
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%   
else
    error('txtltoolbox:txtl_prom_p28_ptet:undefinedmode', ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end 



% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
