% txtl_protein_NRII_H139N.m - protein information for NRII_H139N
% Zoltan A. Tuza Sep 2013

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

function varargout = txtl_protein_NRII_H139N(mode, tube, protein, varargin)

paramObj = txtl_component_config('NRII_H139N');

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode.add_dna_driver, 'Setup Species')
    
    geneData = varargin{1};
    
    defaultBasePairs = {'NRII_H139N','lva','terminator';...
        paramObj.Gene_Length,paramObj.LVA_tag_Length,paramObj.Terminator_Length};
    geneData = txtl_setup_default_basepair_length(tube,geneData,...
        defaultBasePairs);
    
    varargout{1} = geneData;
    
    coreSpecies = {protein.Name};
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');
 
%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif strcmp(mode.add_dna_driver, 'Setup Reactions')
    listOfSpecies = varargin{1};
    p = regexp(listOfSpecies,'^protein NRI(-lva)?$', 'match');
    listOfNRI = vertcat(p{:});
    listOfNRIp = regexprep(listOfNRI, '^protein NRI', 'protein NRI-p');
    for k = 1:size(listOfNRIp,1)
    txtl_addspecies(tube, listOfNRIp{k}, cell(1,1), 'Internal');
    end
    listOfSpecies = get(tube.species,'Name');
    % Need to search again over listOfSPecies because what was just set up may not be all the
    % proteins there are. 
    p = regexp(listOfSpecies,'^protein NRI-p(-lva)?$', 'match');
    listOfProtein = vertcat(p{:});
    replacedProtein = regexprep(listOfProtein, '^protein NRI-p', 'protein NRI');
    
    for k = 1:size(listOfProtein,1)

    % Set up the maturation reaction
    txtl_addreaction(tube,['[' protein.Name '] + ' listOfProtein{k} '  <-> [' protein.Name ':' listOfProtein{k} ']'],...
     'MassAction',{'TXTL_PROT_NRII_H139N_DEPHOSPHORILATION_F',paramObj.generic_rate3;'TXTL_PROT_NRII_H139N_DEPHOSPHORILATION_R',paramObj.generic_rate3});
 
    txtl_addreaction(tube,['[' protein.Name ':' listOfProtein{k} ' ] -> [' protein.Name '] + [' replacedProtein{k} ']'],...
     'MassAction',{'TXTL_PROT_NRII_H139N_DEPHOSPHORILATION',paramObj.generic_rate3});
    end
%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%% 
else
    error('txtltoolbox:txtl_protein_NRII_H139N:undefinedmode', ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end



% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
