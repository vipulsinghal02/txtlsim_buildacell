% txtl_5UTR_utr1b.m - UTR1b ribosome binding site file
% VS, Aug 2017
%
% This file contains a description of the UTR1 ribosome binding site.
% Calling the function txtl_5UTR_UTR1b() will set up the reactions for
% translation with the measured binding rates and translation rates.

% Written by Vipul Singhal 2017
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

function varargout = txtl_5UTR_utr1b(mode, tube, rna, protein, varargin)

    % Create strings for reactants and products
    RNA = ['[' rna.Name ']'];		% RNA species name for reactions
    Ribobound = ['Ribo:' rna.Name];	% Name of bound complex
    
    % importing the corresponding parameters
    paramObj = txtl_component_config('utr1b');

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode.add_dna_driver, 'Setup Species')
    
    utrRbsData = varargin{1};
    defaultBasePairs = {'rbs','spacer';paramObj.RBS_Length,paramObj.spacer_Length};
    utrRbsData = txtl_setup_default_basepair_length(tube,utrRbsData,defaultBasePairs);
    RiboBound = ['Ribo:' rna.Name];
    coreSpecies = {'Ribo',RiboBound};
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');
    varargout{1} = sbioselect(tube, 'Name', RiboBound);
    varargout{2} = utrRbsData;
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%     
elseif strcmp(mode.add_dna_driver, 'Setup Reactions')
    
    % Set up the binding reaction
    txtl_addreaction(tube,['[' rna.Name '] + Ribo <-> [Ribo:' rna.Name ']'],...
        'MassAction',{'TXTL_UTR_UTR1_F',paramObj.Ribosome_Binding_F;
        'TXTL_UTR_UTR1_R',paramObj.Ribosome_Binding_R});
    
    
    
    
    
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%    
else
    error('txtltoolbox:txtl_5UTR_utr1:undefinedmode', ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end    

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:

