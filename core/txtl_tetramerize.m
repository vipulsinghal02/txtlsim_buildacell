% txtl_tetramerize.m - general protein tetramerization
% Zoltan A. Tuza Sep 2012
%
% Reactions for protein tetramerization. 

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

function txtl_tetramerize(mode, tube,protein, varargin)
% function for protein tetramerization.
% tube: sbiomodel object, where the reaction occurs
% protein: SimBiology Species Array
% reacctionRates: 2x1 or 1x2 vector contains the forward and reverse
% reaction rates. (For forward reaction only, set the reverse reaction rete to zero!)
%
% Return: SimBiology Reaction Array

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode.add_dna_driver, 'Setup Species')
    
    coreSpecies = {[protein.Name 'dimer'],[protein.Name 'tetramer']};
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode.add_dna_driver, 'Setup Reactions')
    
%     DONT NEED ALL THIS, The reason is, the parent function is called
%     based on a single protein. it is that very protein we dimerize, and
%     then tetramerize. so dont need to search all proteins. we can assume
%     that the dimer exists. dont need listOfSpecies here. 
%     matchStr = regexp(listOfSpecies,['(^' protein.Name '.*dimer$)'],'tokens','once');
%     listOfproteindimer = vertcat(matchStr{:});
%     if ~isempty(listOfproteindimer)
%         proteintetramerization = true;
%     end
%     
%     if proteintetramerization
%         for i = 1:size(listOfproteindimer,1)
%             proteintetramer = regexprep(listOfproteindimer{i}, 'dimer$', 'tetramer' );
%             Robj = addreaction(tube, ...
%               ['2 [' listOfproteindimer{i} '] <-> [' proteintetramer ']']);
%             Kobj = addkineticlaw(Robj,'MassAction');
%             rN = regexprep(listOfproteindimer{i}, {'( )'}, {''});
%             uniqueNameF = sprintf('TXTL_PROT_TETRAMER_%s_F',rN);
%             uniqueNameR = sprintf('TXTL_PROT_TETRAMER_%s_R',rN);
%             set(Kobj, 'ParameterVariableNames', {uniqueNameF, uniqueNameR});
%         end
%     end

 reactionRate = varargin{1};
    
    if isempty(reactionRate)
        error('txtltoolbox:txtl_tetramerize:unspecifiedRR', ...
        'Please specify default tetramerization reaction rates as an input vector');
    else
       Robj = addreaction(tube, ['2 [' protein.Name 'dimer] <-> [' protein.Name 'tetramer]']); 
       Kobj = addkineticlaw(Robj,'MassAction');
       rN = regexprep(protein.Name, {'( )'}, {''});
       uniqueNameF = sprintf('TXTL_PROT_TETRAMER_%s_F',rN);
       uniqueNameR = sprintf('TXTL_PROT_TETRAMER_%s_R',rN);
       Pobjf = addparameter(Kobj, uniqueNameF, reactionRate(1)); 
       Pobjr = addparameter(Kobj, uniqueNameR, reactionRate(2));
       set(Kobj, 'ParameterVariableNames', {uniqueNameF, uniqueNameR});
    end
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%    
else
    error('txtltoolbox:txtl_tetramerize:undefinedmode', ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end   

end



