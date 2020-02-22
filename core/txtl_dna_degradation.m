% txtl_dna_degradation.m - general protein degradation model
% Vipul Singhal Sep 2012
%
% This file contains a description of the dna degradation reactions.
% Calling the function txtl_dna_degradation will set up the reactions for
% the degradation of any linear dna by RecBCD.
% reactionRates is a 3x1 or 1x3 vector. 

%
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


function txtl_dna_degradation(mode, tube,dna,varargin)
% function for dna degradation.
% tube: sbiomodel object, where the reaction occurs
% dna: simBiology species object
% reacctionRate: degradation rate
%
%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Species %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode.add_dna_driver, 'Setup Species')
    
    coreSpecies = {'RecBCD',[dna.Name ':RecBCD']};
    % empty cellarray for amount => zero amount
    txtl_addspecies(tube, coreSpecies, cell(1,size(coreSpecies,2)), 'Internal');
    
%%%%%%%%%%%%%%%%%%% DRIVER MODE: Setup Reactions %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode.add_dna_driver,'Setup Reactions')
    
    reactionRates = varargin{1};    
    
    switch length(reactionRates)
        case 3
           Robj1 = addreaction(tube, [ dna.Name ' + RecBCD <-> [' dna.Name ':RecBCD]']);
           Kobj1 = addkineticlaw(Robj1,'MassAction');
           Pobj1 = addparameter(Kobj1,  'kf_complex', reactionRates(1));
           Pobj1r = addparameter(Kobj1, 'kr_complex', reactionRates(2));
           set(Kobj1, 'ParameterVariableNames', {'kf_complex', 'kr_complex'});

           Robj2 = addreaction(tube, ['[' dna.Name ':RecBCD] -> RecBCD']);
           Kobj2 = addkineticlaw(Robj2,'MassAction');
           Pobj2 = addparameter(Kobj2,  'kf_deg', reactionRates(3));
           set(Kobj2, 'ParameterVariableNames', 'kf_deg'); 
        case 2
           Robj1 = addreaction(tube, [ dna.Name ' + RecBCD -> [' dna.Name ':RecBCD]']);
           Kobj1 = addkineticlaw(Robj1,'MassAction');
           Pobj1 = addparameter(Kobj1,  'kf_complex', reactionRates(1));
           set(Kobj1, 'ParameterVariableNames', {'kf_complex', 'kr_complex'});

           Robj2 = addreaction(tube, ['[' dna.Name ':RecBCD] -> RecBCD']);
           Kobj2 = addkineticlaw(Robj2,'MassAction');
           Pobj2 = addparameter(Kobj2,  'kf_deg', reactionRates(2));
           set(Kobj2, 'ParameterVariableNames', 'kf_deg');    
        case 1
           Robj1 = addreaction(tube, [ dna.Name ' + RecBCD -> RecBCD']);
           Kobj1 = addkineticlaw(Robj1,'MassAction');
           Pobj1 = addparameter(Kobj1,  'kf_deg', reactionRates(1));
           set(Kobj1, 'ParameterVariableNames', 'kr_deg');
    end
%%%%%%%%%%%%%%%%%%% DRIVER MODE: error handling %%%%%%%%%%%%%%%%%%%%%%%%%%%    
else
    error('txtltoolbox:txtl_dna_degradation:undefinedmode', ...
      'The possible modes are ''Setup Species'' and ''Setup Reactions''.');
end 

end
