% txtl_setup_default_basepair_length
% ZT, Dec 2012
%
% In case the user won't provide any promoter size this function sets it up
% based on the defaultbasepair array.
% Every promoter and utr_rbs files should call this file
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


function ObjectData = txtl_setup_default_basepair_length(tube,ObjectData,defaultBasePairs)


    %TODO! 12/08/12 zoltuz - later on deafaultBasePairs should come from
    %TODO the config files - that's way tube is provided 
    noLength = cellfun('isempty', ObjectData(2,:));
    elementInd = find(noLength>0); % indexes where is no length
    % find the elements with missing length
    [~,Aind,Bind] = intersect(defaultBasePairs(1,:),ObjectData(1,elementInd)); 
    
    % multiple lhs - multiple rhs -> no way to avoid for cycle
    if ~isempty(Aind)
        for k = 1:size(Aind,2)
            ObjectData{2,elementInd(Bind(k))} = defaultBasePairs{2,Aind(k)};
        end
    end



end