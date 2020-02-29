% txtl_prom_OR1OR2.m - promoter information for pOR1OR2 promoter
% VS July 2014
%
% This file contains a description of the ptet promoter.
% Calling the function txtl_prom_pOR1OR2() will set up the reactions for
% transcription with the measured binding rates and transription rates.

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

function varargout = txtl_prom_pOR2OR1(mode, tube, dna, rna,varargin)

% Create strings for reactants and products
DNA = ['[' dna.Name ']'];		% DNA species name for reactions
RNA = ['[' rna.Name ']'];		% RNA species name for reactions
RNAP = 'RNAP';			% RNA polymerase name for reactions
RNAPbound = ['RNAP:' dna.Name];
% importing the corresponding parameters
paramObj = txtl_component_config('pOR2OR1');


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
    defaultBasePairs = {'pOR2OR1','junk','thio';...
        paramObj.Promoter_Length,paramObj.Junk_Length,paramObj.Thio_Length};
    promoterData = txtl_setup_default_basepair_length(tube,promoterData,...
        defaultBasePairs);
    
    varargout{1} = promoterData;
    
    coreSpecies = {RNAP,RNAPbound};
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');
%     if mode.utr_attenuator_flag
%         txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, RNAPbound, prom_spec, rbs_spec, gene_spec );
%     else
        txtl_transcription(mode, tube, dna, rna, RNAP, RNAPbound);
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
    % Parameters that describe this promoter
    parameters = {'TXTL_POR2OR1_RNAPbound_F',paramObj.RNAPbound_Forward;...
        'TXTL_POR2OR1_RNAPbound_R',paramObj.RNAPbound_Reverse};
    % Set up binding reaction
    txtl_addreaction(tube,[DNA ' + ' RNAP ' <-> [' RNAPbound ']'],...
        'MassAction',parameters);
   
%     if mode.utr_attenuator_flag
%         txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, RNAPbound, prom_spec, rbs_spec, gene_spec );
%     else
        txtl_transcription(mode, tube, dna, rna, RNAP, RNAPbound);
%     end
    
    
    matchStr = regexp(listOfSpecies,'(^protein cLambda.*dimer$)','tokens','once');
    listOfcLambdadimer = vertcat(matchStr{:});
    
    
    % repression of ptet by tetR dimer
    if ~isempty(listOfcLambdadimer)
        for k = 1:size(listOfcLambdadimer,1)
            txtl_addreaction(tube,...
                [DNA ' + ' listOfcLambdadimer{k} ' <-> [' dna.name ':' listOfcLambdadimer{k} ']'],...
                'MassAction',{'pOR2OR1_sequestration_F',getDNASequestrationRates(paramObj,'F');...
                'pOR1OR2_sequestration_R',getDNASequestrationRates(paramObj,'R')});
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error('txtltoolbox:txtl_prom_pOR2OR1:undefinedmode', ...
        'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end
end


% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
