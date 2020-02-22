% txtl_prom_plas.m - promoter information for plas
% VS Dec 2013
%
% This file contains a description of the plasR promoter, which is activated by the lasR protein. .
% Calling the function txtl_prom_plas() will set up the reactions for
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

function varargout= txtl_prom_plas(mode, tube, dna, rna, varargin)


% Create strings for reactants and products
DNA = ['[' dna.Name ']'];		% DNA species name for reactions
RNA = ['[' rna.Name ']'];		% RNA species name for reactions
RNAP = 'RNAP';			% RNA polymerase name for reactions
RNAPbound = ['RNAP:' dna.Name];

% importing the corresponding parameters
paramObj = txtl_component_config('plas');

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
    defaultBasePairs = {'plas','junk','thio';150,500,0};
    promoterData = txtl_setup_default_basepair_length(tube,promoterData,...
        defaultBasePairs);
    
    varargout{1} = promoterData;
    
    coreSpecies = {RNAP,RNAPbound};
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');
    %{
    % !TODO: VS Aug 2014: dont yet know what form of P3 exists: is lasR lva/ssrA bound or not?
    % !TODO: SOLUTION VS Aug 2017: have a multipass method for generating
    the toolbox equations. basically keep generating species and equations
    until no more changes occur. 
    
    % VS Aug 2014: Comment txtl_transcription out for now. the species and
    reactions all get set up in the setup reactions phase, when we have
    sufficient information about which versionof the protein will be
    present. later we will generalize the toolbox to be a multipass method
    which never has any of these probelms. 
    
    
    %     if mode.utr_attenuator_flag
    % %         txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP,RNAPbound, prom_spec, rbs_spec, gene_spec );
    %         txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' P3 ],prom_spec, rbs_spec, gene_spec,{P3} );
    %     else
    % %         txtl_transcription(mode, tube, dna, rna, RNAP,RNAPbound);
    %         txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' P3 ],{P3});
    %     end
    %}
    
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
    
    %   leaky expression (remove and make conditional on the leaky
    %   expression mode?)
        parameters = {'TXTL_PLAS_RNAPbound_F',paramObj.RNAPbound_Forward;...
                      'TXTL_PLAS_RNAPbound_R',paramObj.RNAPbound_Reverse};
        txtl_addreaction(tube,[DNA ' + ' RNAP ' <-> [' RNAPbound ']'],...
                'MassAction',parameters);
    
    p = regexp(listOfSpecies,'^OC12HSL:protein lasR(-lva)?$', 'match');
    TF = vertcat(p{:});
    for k = 1:size(TF,1)
        setup_promoter_reactions(mode, tube, dna, rna, RNAP,RNAPbound,...
            prom_spec, rbs_spec, gene_spec,TF{k}, paramObj)
    end
    
    %%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error('txtltoolbox:txtl_prom_plas:undefinedmode', ...
        'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end
end

function setup_promoter_reactions(mode, tube, dna,rna, RNAP, RNAPbound,...
    prom_spec, rbs_spec, gene_spec, TF, paramObj)

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
    'TXTL_PLAS_TFBIND_R',paramObj.activation_R});

% 
% if mode.utr_attenuator_flag
%     mode.add_dna_driver = 'Setup Species';
%     txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' TF],prom_spec, rbs_spec, gene_spec,{TF} );
%     mode.add_dna_driver = 'Setup Reactions';
%     txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' TF],prom_spec, rbs_spec, gene_spec,{TF} );
% else
    mode.add_dna_driver = 'Setup Species';
    txtl_transcription(mode, tube, dna, rna, RNAP,RNAPbound);
    txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' TF ],{TF});
    mode.add_dna_driver = 'Setup Reactions';
    txtl_transcription(mode, tube, dna, rna, RNAP,RNAPbound);
    txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' TF ],{TF});
% end

end



% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
