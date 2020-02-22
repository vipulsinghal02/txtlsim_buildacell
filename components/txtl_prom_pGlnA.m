% txtl_prom_pGlnA.m - promoter information for glnA promoter file
% Zoltan Tuza, Sep 2013
%
% Written by Zoltan Tuza, Sep 2013
%
% Copyright (c) 2013 by California Institute of Technology
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

function varargout= txtl_prom_pGlnA(mode, tube, dna, rna, varargin)
% Create strings for reactants and products
DNA = ['[' dna.Name ']'];		% DNA species name for reactions
RNA = ['[' rna.Name ']'];		% RNA species name for reactions
RNAP = 'RNAP54';			% RNA polymerase name for reactions
RNAPbound = ['RNAP54:' dna.Name];
P1 = 'protein sigma54';
% importing the corresponding parameters
paramObj = txtl_component_config('pGlnA');

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
    defaultBasePairs = {'pGlnA','junk','thio';...
        paramObj.Promoter_Length,paramObj.Junk_Length,paramObj.Thio_Length};
    promoterData = txtl_setup_default_basepair_length(tube,promoterData,...
        defaultBasePairs);
    
    varargout{1} = promoterData;
    
    coreSpecies = {RNAP,RNAPbound,P1};
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');
    
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
    parameters = {'TXTL_pGlnA_RNAPbound_F',paramObj.RNAPbound_Forward;...
        'TXTL_pGlnA_RNAPbound_R',paramObj.RNAPbound_Reverse};
    txtl_addreaction(tube,[DNA ' + ' RNAP ' <-> [' RNAPbound ']'],...
        'MassAction',parameters);
    
    p = cellfun(@(x) regexp(x,'^protein NRI-p(-lva)?$', 'match'), listOfSpecies, 'UniformOutput', false);
    listOfProtein = vertcat(p{:});
    
    % NRI-p
    for k = 1:size(listOfProtein,1)
        txtl_addreaction(tube,[dna.Name ' + ' listOfProtein{k} ' <-> ' dna.Name ':' listOfProtein{k} ],...
            'MassAction',{'TXTL_DNA_NRI-p_F',paramObj.activation_F;'TXTL_DNA_NRI-p_R',paramObj.activation_R});
        txtl_addreaction(tube,[RNAPbound  ' + ' listOfProtein{k} ' <-> [' RNAPbound ':' listOfProtein{k} ']'],...
            'MassAction',{'TXTL_RNAPbound_NRI-p_F',paramObj.activation_F;'TXTL_RNAPbound_NRI-p_R',paramObj.activation_R});
        txtl_addreaction(tube,[dna.Name ':' listOfProtein{k} ' + ' RNAP ' <-> [' RNAPbound ':' listOfProtein{k} ']' ],...
            'MassAction',{'TXTL_DNA_NRI-p_RNAPbound_F',paramObj.RNAPbound_Forward_actv;'TXTL_DNA_NRI-p_RNAPbound_R',paramObj.RNAPbound_Reverse_actv});
        
        
%         if mode.utr_attenuator_flag
%             mode.add_dna_driver = 'Setup Species';
%             txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP,RNAPbound, prom_spec, rbs_spec, gene_spec );
%             txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' listOfProtein{k} ],prom_spec, rbs_spec, gene_spec,{listOfProtein{k} } );
%             mode.add_dna_driver = 'Setup Reactions';
%             txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP,RNAPbound, prom_spec, rbs_spec, gene_spec );
%             txtl_transcription_RNAcircuits(mode, tube, dna, rna, RNAP, [RNAPbound ':' listOfProtein{k} ],prom_spec, rbs_spec, gene_spec,{listOfProtein{k}} );
%         else
            mode.add_dna_driver = 'Setup Species';
            txtl_transcription(mode, tube, dna, rna, RNAP,RNAPbound);
            txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' listOfProtein{k} ],{listOfProtein{k}});
            mode.add_dna_driver = 'Setup Reactions';
            txtl_transcription(mode, tube, dna, rna, RNAP,RNAPbound);
            txtl_transcription(mode, tube, dna, rna, RNAP,[RNAPbound ':' listOfProtein{k} ],{listOfProtein{k}});
%         end
    end
 
    
    %%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error('txtltoolbox:txtl_prom_pGlnA:undefinedmode', ...
        'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end



% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
