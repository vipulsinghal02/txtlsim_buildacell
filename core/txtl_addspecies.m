function simBioSpecies = txtl_addspecies(tube, name, amount, varargin)
%TXTL_ADDSPECIES   Add one or more molecular species to a tube
%
% Sobj = TXTL_ADDSPECIES(tube, name, amount) adds a molecule to a
% tube, in the gen amount (in nM).  The species can be a new species or
% one that already exists (in which case its concentration is added to
% what is already present).

% Written by Richard Murray, 11 Sep 2012
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

add_dna_mode = struct('add_dna_driver', {'Setup Species'});


% checks
if (length(amount) > 1)
    assert(length(name) == length(amount));
end

if ~iscell(amount)
    amount = num2cell(amount);
end

if ischar(name)
    name = {name};
end

%%%%% Handling multiple compartments 
if size(tube.Compartments,1) > 1
    % test for valid compartment name in species string
    [matchstr splitstr] = regexp(name,'\.','split');
    if sum(find(cellfun(@isempty,splitstr) == 1)) > 0
        error('Model has multiple compartmerts, but no a valid compartment ID was given!');
    end
    compList = get(tube.Compartments,'name');
    
    for k = 1:length(matchstr);
        matchRes = findStringInAList(compList,matchstr{k}{1});
        name{k} = matchstr{k}{2};
        ind = find(matchRes > 0);
        if ~isempty(ind)
            compInd(k) = ind;
        else
            compInd(k) = 0;
        end
    end
    
    if sum(find(compInd == 0)) ~= 0
        error('Model has multiple compartmerts, but no (not a) valid compartment ID was given!');
    else
        currentComp = compList(compInd);
        % addspecies to the corresponding compartments
        for k=1:size(compInd,1)
            index = findspecies(tube.Compartments(compInd(k)), name{k});
            index = num2cell(index);
            addOneSpecie(tube.Compartments(compInd(k)),name{k},amount{k},index{k});
        end
        % TODO zoltuz 30/7/13 possible error with multi compartment species
        % e.g. NTP,IPTG exist in multiple compartments
        indexPost = findspecies(tube, name);
        simBioSpecies = tube.Species(indexPost);
        return
    end
else
    index = findspecies(tube, name);
end


index = num2cell(index);

cellfun(@(x,y,z) addOneSpecie(tube,x,y,z,'no_output'),name,amount,index,'UniformOutput',false);
indexPost = findspecies(tube, name);
simBioSpecies = tube.Species(indexPost);

% has protein been added?
 if isempty(varargin) && sum(cellfun(@(x) isempty(strfind(x,'protein')),name)) >0
       txtl_setup_new_protein_added(tube,add_dna_mode);
 end
end


function varargout = addOneSpecie(tube,name,amount,index,varargin)
% if amount wasn't specified then make it zero
if isempty(amount)
    amount = 0;
end

if (index == 0)
    varargout{1} = addspecies(tube, name, amount);
else
    tube.Species(index).InitialAmount = ...
        tube.Species(index).InitialAmount + amount;
    tube.Species(index);
    varargout{1} = tube.Species(index);
end

end


% Automatically use MATLAB mode in Emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
