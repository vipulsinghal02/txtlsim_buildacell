% txtl_dimerize.m - general protein dimerization
% Zoltan A. Tuza Sep 2012
%
%
% This file contains a description of the protein produced by tetR.
% Calling the function txtl_protein_tetR() will set up the reactions for
% sequestration by the inducer aTc.

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

function txtl_dimerize(mode, tube,protein,varargin)
% function for protein dimerization.
% tube: sbiomodel object, where the reaction occurs
% protein: SimBiology Species Array
% reacctionRates: 2x1 or 1x2 vector contains the forward and reverse
% reaction rates. (For forward reaction only, set the reverse reaction rete to zero!)
%

%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode.add_dna_driver, 'Setup Species')
    
    txtl_addspecies(tube, [protein.Name 'dimer'], 0, 'Internal');
      
%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%       
elseif strcmp(mode.add_dna_driver, 'Setup Reactions')
    
    reactionRate = varargin{1};
    regIX = regexp(protein.Name,'^protein ', 'end');
    protstr = protein.Name((regIX+1):end);
    if isempty(reactionRate)
        error('txtltoolbox:txtl_dimerize:unspecifiedRR', ...
            'Please specify dimerization reaction rates as an input vector');
    elseif reactionRate(2) == 0 || length(reactionRate) == 1
       Robj = addreaction(tube, ['2 [' protein.Name '] -> [' protein.Name 'dimer]']);
       Kobj = addkineticlaw(Robj,'MassAction');
       Pobj = addparameter(Kobj,  ['TXTL_DIMER_' protstr '_' 'F'], reactionRate(1));
       set(Kobj, 'ParameterVariableNames',['TXTL_DIMER_' protstr '_' 'F']);
    else
       Robj = addreaction(tube, ['2 [' protein.Name '] <-> [' protein.Name 'dimer]']); 
       Kobj = addkineticlaw(Robj,'MassAction');
       Pobjf = addparameter(Kobj, ['TXTL_DIMER_' protstr '_' 'F'], reactionRate(1));
       Pobjr = addparameter(Kobj, ['TXTL_DIMER_' protstr '_' 'R'], reactionRate(2));
       set(Kobj, 'ParameterVariableNames', {['TXTL_DIMER_' protstr '_' 'F'],...
           ['TXTL_DIMER_' protstr '_' 'R']});
    end
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error('txtltoolbox:txtl_dimerize:undefinedmode', ...
       'The possible modes are ''Setup Species'' and ''Setup Reactions''');
end


end



