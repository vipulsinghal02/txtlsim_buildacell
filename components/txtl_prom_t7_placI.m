% txtl_prom_t7_placI.m - promoter information for t7 and placI combinatorial promoter
% Zoltan Tuza, July 2013
%
% This file contains a description of the t7 and placI combinatorial promoter.
% Calling the function txtl_prom_t7_placI() will set up the reactions for
% transcription with the measured binding rates and transription rates.
% 
% 

% Written by Zoltan Tuza, July 2013
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


%% Promoter activity
%
% T7RNAP + DNA -> T7RNAP:DNA:NTP -> T7RNAP + DNA + mRNA (High)
% RNAP70 + DNA -> RNAP70:DNA:NTP -> RNAP70 + DNA + mRNA (low)
% T7RNAP + DNA:protein LacI -> T7RNAP:DNA:protein LacI -> T7RNAP +
% DNA:protein LacI +  mRNA (leaky) 
%

function varargout= txtl_prom_t7_placI(mode, tube, dna, rna, varargin)


    % Create strings for reactants and products
    DNA = ['[' dna.Name ']'];		% DNA species name for reactions
    RNA = ['[' rna.Name ']'];		% RNA species name for reactions
    RNAP = 'RNAP70';			% RNA polymerase name for reactions
    RNAPbound = ['RNAP70:' dna.Name];
    T7RNAP = 'protein t7RNAP';			% RNA polymerase name for reactions
    T7RNAPbound = ['protein t7RNAP:' dna.Name];	% Name of bound complex
    P2 = 'protein LacI';
    
    % importing the corresponding parameters
    paramObj = txtl_component_config('t7_lacI');
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode.add_dna_driver, 'Setup Species')

    
    promoterData = varargin{1};
    defaultBasePairs = {'t7_lacI','junk','thio';150,500,0};
    promoterData = txtl_setup_default_basepair_length(tube,promoterData,...
        defaultBasePairs);
    
    varargout{1} = promoterData;
    
    coreSpecies = {RNAP,RNAPbound,T7RNAP,T7RNAPbound,P2};
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)),'Internal');
    
    % nominal transcription
    txtl_transcription(mode, tube, dna, rna, RNAP,RNAPbound);
     % T7RNAP transcription
    txtl_transcription(mode, tube, dna, rna, T7RNAP,T7RNAPbound);

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif strcmp(mode.add_dna_driver,'Setup Reactions')
    
    % Parameters that describe this promoter
    parameters = {'TXTL_RNAPbound_F',paramObj.RNAPbound_Forward;...
                  'TXTL_RNAPbound_R',paramObj.RNAPbound_Reverse};
    % Set up binding reaction
    txtl_addreaction(tube,[DNA ' + ' RNAP ' <-> [' RNAPbound ']'],...
        'MassAction',parameters);
    %
    % nominal transcription
    txtl_transcription(mode, tube, dna, rna, RNAP, RNAPbound);
    % T7RNAP transcription
    txtl_transcription(mode, tube, dna, rna, T7RNAP,T7RNAPbound);

    % 
    %
    Robj3 = addreaction(tube, [dna.Name ' + ' P2 ' <-> ' dna.Name ':' P2 ]);
    Kobj3 = addkineticlaw(Robj3, 'MassAction');
    Pobj3f = addparameter(Kobj3, 'kf', 2.86e-3);
    Pobj3r = addparameter(Kobj3, 'kr', 0.11e-4);
    set(Kobj3, 'ParameterVariableNames', {'kf', 'kr'});
    % 
%     Robj4 = addreaction(tube, [RNAPbound  ' + ' P2 ' <-> ' RNAPbound ':' P2 ]);
%     Kobj4 = addkineticlaw(Robj4, 'MassAction');
%     Pobj4f = addparameter(Kobj4, 'kf', 2.86e-3);
%     Pobj4r = addparameter(Kobj4, 'kr', 0.11e-4);
%     set(Kobj4, 'ParameterVariableNames', {'kf', 'kr'});
    % 
    % Set up binding reaction for LacI
    Robj2 = addreaction(tube, [dna.Name ':' P2 ' + ' T7RNAP ' <-> [' T7RNAPbound ':' P2 ']' ]);
    Kobj2 = addkineticlaw(Robj2, 'MassAction');
    Pobj2f = addparameter(Kobj2, 'kf', paramObj.RNAPbound_Forward);
    Pobj2r = addparameter(Kobj2, 'kr', paramObj.RNAPbound_Reverse*1000);
    set(Kobj2, 'ParameterVariableNames', {'kf', 'kr'});

    txtl_transcription(mode, tube, dna, rna, T7RNAP,[T7RNAPbound ':' P2 ],{P2});
    
    
    matchStr = regexp(listOfSpecies,'(^protein lacI.*tetramer$)','tokens','once'); % ^ matches RNA if it occust at the beginning of an input string
    listOflacItetramers = vertcat(matchStr{:});
    % repression of placI by lacItetramer
    if ~isempty(listOflacItetramers)
        for i = 1:size(listOflacItetramers,1)
             txtl_addreaction(tube,...
                [DNA ' + ' listOflacItetramers{i} ' <-> [' dna.name ':' listOflacItetramers{i} ']'],...
            'MassAction',{'placI_sequestration_F',getDNASequestrationRates(paramObj,'F');...
                          'placI_sequestration_R',getDNASequestrationRates(paramObj,'R')});
        end
    end
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%   
else
    error('txtltoolbox:txtl_prom_t7_lacI:undefinedmode', ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end 



% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
