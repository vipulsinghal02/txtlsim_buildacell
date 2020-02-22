% txtl_protein_lacI3.m - protein information for lacI
%
% This file contains a description of the protein produced by tetR.
% Calling the function txtl_protein_tetR() will set up the reactions for
% sequestration by the inducer aTc.

% Sep 2012
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

function varargout = txtl_protein_lacI3(mode, tube, protein, varargin)

% importing the corresponding parameters
paramObj = txtl_component_config('lacItoggleSwitch2'); 

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode.add_dna_driver, 'Setup Species')

    geneData = varargin{1};
    defaultBasePairs = {'lacI3','lva','terminator';...
        paramObj.Gene_Length,paramObj.LVA_tag_Length,paramObj.Terminator_Length};
    geneData = txtl_setup_default_basepair_length(tube,geneData,...
        defaultBasePairs);
    
    varargout{1} = geneData;
    
    coreSpecies = {'IPTG'}; %,['IPTG:' protein.Name] 
    % !TODO: VS 3 April 13. Not sure, but i think IPTG does NOT bind with the multimerized protein. 
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');
    
    
    % call additional functions to setup any other relevant species (like
    % multimers)
    txtl_dimerize(mode, tube, protein); 
    txtl_tetramerize(mode, tube, protein);
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode.add_dna_driver, 'Setup Reactions')
   
    [~,listOfSpecies] = getstoichmatrix(tube);
    
    % Set up the binding reaction for all protein variants
    matchStr = regexp(listOfSpecies,'(^protein lacI3.*tetramer$)','tokens','once'); 
    listOftetramers = vertcat(matchStr{:});

    for k = 1:size(listOftetramers,1)
        txtl_addreaction(tube, ...
         ['[' listOftetramers{k} '] + IPTG <-> [IPTG:' listOftetramers{k} ']'],...
         'MassAction',{'TXTL_INDUCER_LACI_IPTG_F',paramObj.Protein_Inducer_Forward;...
                       'TXTL_INDUCER_LACI_IPTG_R',paramObj.Protein_Inducer_Reverse});
    end
    
     % degrade the IPTG inducer
    txtl_addreaction(tube,'IPTG -> null',...
     'MassAction',{'TXTL_INDUCER_DEGRADATION_IPTG',paramObj.Inducer_Degradation});
    

    % set up a reaction for protein dimerization
    txtl_dimerize(mode, tube,protein, ...
        [paramObj.Dimmerization_Forward, paramObj.Dimmerization_Reverse]);

    %Set up tetramerization
    txtl_tetramerize(mode, tube,protein,...
        [paramObj.Tetramerization_Forward, paramObj.Tetramerization_Reverse]);

%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%    
else
    error('txtltoolbox:txtl_protein_lacI3:undefinedmode', ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end    





% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
