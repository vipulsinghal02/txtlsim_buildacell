% txtl_prom_ptet2.m - promoter information for ptet promoter
% VS Sep 2012
%
%! TODO: Header needs to be updated (RMM didn't write this file; 29 Sep 2012)
%! TODO: This file is not properly named; promoter is ptet, not ptet2 ??
% 
% This file contains a description of the ptet promoter.
% Calling the function txtl_prom_ptet2() will set up the reactions for
% transcription with the measured binding rates and transription rates.
% The binding of the promoter to the lacI repressor is used in the
% gen_switch example. (the original file, txtl_prom_ptet.m is for negative
% autoregulation, where tetR represses itself. 

% Written by Vipul Singhal, Sep 2012
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

function Rlist = txtl_prom_ptet2(tube, dna, rna)

% Parameters that describe this promoter
%! TODO: replace these values with correct values
kf_ptet = log(2)/0.1;			% 100 ms bind rate
kr_ptet = 10 * kf_ptet;			% Km of 10 (same as p70, from VN)
ktx_ptet = log(2)/(rna.UserData/30);	% 30 base/second transcription

% Create strings for reactants and products
DNA = ['[' dna.Name ']'];		% DNA species name for reactions
RNA = ['[' rna.Name ']'];		% RNA species name for reactions
RNAP = 'RNAP70';			% RNA polymerase name for reactions
RNAPbound = ['RNAP70:' dna.Name];

% Set up binding reaction
Robj1 = addreaction(tube, [DNA ' + ' RNAP ' <-> [' RNAPbound ']']);
Kobj1 = addkineticlaw(Robj1, 'MassAction');
Pobj1f = addparameter(Kobj1, 'kf', kf_ptet);
Pobj1r = addparameter(Kobj1, 'kr', kr_ptet);
set(Kobj1, 'ParameterVariableNames', {'kf', 'kr'});

%
% Now put in the reactions for the utilization of NTPs
% Use an enzymatic reaction to proper rate limiting
%

txtl_transcription(tube, dna, rna, RNAP, RNAPbound);

%
% Add reactions for sequestration of promoter by tetRdimer 
%

%! TODO: Looks like lacI is used instead of tetR??? (RMM, 29 Sep 2012)
%! TODO: I don't see the dimer binding; just the monomer??? (RMM, 29 Sep 2012)
% VS, 1 Oct 2012: About ptet2: I need to specify the reaction between the promoter ptet and 
% its repression by lacI. in the original ptet file, tetR represses ptet, 
% but for the purposes of the genetic switch, we need the lacI dimer to 
% repress the ptet promoter. 

% So depending on the application (negautoreg or gen_switch), we may want 
% the promoter to have different interactions with proteins. A possible 
% solution is that we specify these interactions in a separate file, and not 
% in the promoter file. 

% I was not aware that dimerization is a requirement for the repression. 
% In any case, the monomer was temporary, because its concentration was high 
% enough for the switching due to the repression to be effective. We can replace 
% this with the dimer, but that will require higher dimerization rates, and 
% I wanted to discuss how to obtain some of these reaction rates before changing that. 

kf_lacI = 4; kr_lacI = 0.1;		% 
Robj4 = addreaction(tube, ...
  [DNA ' + [protein lacI] <-> [DNA tetR:protein lacI]']);
Kobj4 = addkineticlaw(Robj4,'MassAction');
Pobj4 = addparameter(Kobj4, 'k4', kf_lacI);
Pobj4r = addparameter(Kobj4, 'k4r', kr_lacI);
set(Kobj4, 'ParameterVariableNames', {'k4', 'k4r'});

Rlist = [Robj1, Robj4];

% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
