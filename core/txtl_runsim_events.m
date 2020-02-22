% txtl_runsim_events.m - template file for MATLAB functions
% Vipul Singhal Nov 2012
%
% This file is a template that includes the boilerplate for creating a 
% MATLAB function with all of the BSD licensing information at the top.


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


function [t_ode, x_ode, simData] = txtl_runsim_events(varargin)
        
Mobj = varargin{1};
configsetObj = varargin{2};
eventTriggers = varargin{3};
eventFcns = varargin{4};
evt = cell(1, length(eventTriggers));
for i = 1:length(eventTriggers)
evt{i} = addevent(Mobj, eventTriggers{i}(1:end), eventFcns{i}(1:end));
end

if isempty(configsetObj)
    configsetObj = getconfigset(Mobj);
end

set(configsetObj.RuntimeOptions, 'StatesToLog', 'all');

simData = sbiosimulate(Mobj);

[t_ode,x_ode] = deal(simData.Time, simData.Data);

end



% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
