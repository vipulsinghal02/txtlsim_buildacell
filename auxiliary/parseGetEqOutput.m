% this functions parse the Simbiology's getequations out. it generates an
% ode file and parameter, initial value data strcutrure. 
% Input: Simbiology Model
% Output: data structure with paramteres, initial values and a path to the
% ode file


% Written by Zoltan A Tuza, Aug 2013
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
%%
function dataOut = parseGetEqOutput(Mobj) 


if ismethod(Mobj,'getequations')
    str = getequations(Mobj);
else 
    error('getequations function does not exist in the currently installed Simbiology version, upgrade your matlab version at least to Matalb 2012b!');
end


tic 
linebyline = textscan(str, '%s', 'delimiter', sprintf('\n'));

% find categories (Fluxes,Parameter Values,Initial Conditions,etc)
div = regexp(linebyline{1},'^(.*):$','tokens','once');
% save the line number where a category is located
divPos = find(~cellfun('isempty',div) >0);

% most of the string output is expression type data, but some are table type
tableType = {'Events'};

for k = 1:size(divPos,1)
    processedOutput{1,k} = div{divPos(k)};
    
    % cut out the right range between categories
    if k < size(divPos,1)
        range = divPos(k)+1:divPos(k+1)-2;
    else
        range = divPos(k)+1:size(linebyline{1},1);
    end
    
    % is the current category a 'expression' type? check table type list
    if strcmp(processedOutput{1,k},tableType) == 0
        tmp = parseExpression(linebyline{1}(range),' = ');
    else
        tmp = parseTable(linebyline{1}(range),'=');
    end
    processedOutput{2,k} = vertcat(tmp{:});
    
end


%% write ode file header
fileName = Mobj.Name;
filePath = ['tmp/' fileName '.m'];
fid = fopen(filePath,'w');

fprintf(fid,'function dx = %s(time,x,p)\n\n\n',fileName);

% find fluxes
fluxesCell = cellfun(@(x) strcmp(x,'Fluxes'),processedOutput(1,:));



% process ode equations
odeEqcell = cellfun(@(x) strcmp(x,'ODEs'),processedOutput(1,:));
numOfSpecies = size(processedOutput{2,odeEqcell},1);
processedOutput{2,odeEqcell}(:,3) = cellfun(@(z) sprintf('x(%d)',z),num2cell(1:numOfSpecies),'UniformOutput',false);
r = regexp(processedOutput{2,odeEqcell}(:,1),'^d\((.*)\)/dt','tokens','once');
processedOutput{2,odeEqcell}(:,1) = vertcat(r{:});
% process parameters
paramCell = cellfun(@(x) strcmp(x,'Parameter Values'),processedOutput(1,:));
numOfParameters = size(processedOutput{2,paramCell},1);
processedOutput{2,paramCell}(:,3) = cellfun(@(z) sprintf('p(%d)',z),num2cell(1:numOfParameters),'UniformOutput',false);
% find events
eventsCell = cellfun(@(x) strcmp(x,'Events'),processedOutput(1,:));
processedOutput{2,eventsCell}(:,2) = regexp(processedOutput{2,eventsCell}(:,2),' = ','split','once');
% find events in the parameter vector, insert if necessary
for k=1:size(processedOutput{2,eventsCell}(:,2),1)
    oo = processedOutput{2,eventsCell}(:,2);
    idx = find(strcmp(oo{k,1}{1},processedOutput{2,paramCell}(:,1))>0);
    if isempty(idx)
        % add the new parameter and its value along with the incremented
        % identifier 
        pID = sprintf('p(%d)',size(processedOutput{2,paramCell}(:,1),1)+1);
        processedOutput{2,paramCell}(end+1,:) = {oo{k,1}{1},oo{k,1}{2},pID};
        oo{k,1}{1} = pID;
        processedOutput{2,eventsCell}(:,2) = oo;
    else
        oo{k,1}{1} = processedOutput{2,paramCell}{idx,3};
        processedOutput{2,eventsCell}(:,2) = oo;
    end
end

%%
% multi string replace

for k = 1:size(processedOutput{2,fluxesCell}(:,1))
    currentFlux = processedOutput{2,fluxesCell}(k,2);
    
    
    [matchstr,splitstr] = regexp(currentFlux,'\*|\+|(?<=\w|])\-(?=\w|[)|\/|\(|\)','match','split');
    
    splitstr = splitstr{1};
    
    % replacing state variables
    oo = cellfun(@(y) find(strcmp(y,processedOutput{2,odeEqcell}(:,1))>0),splitstr,'UniformOutput',false);
    ind= ~cellfun('isempty',oo);
    elementPos = find(ind >0);
    sel = cell2mat(oo(ind));
    splitstr(elementPos) = processedOutput{2,odeEqcell}(sel,3);
    
    % replacing parameters 
    oo = cellfun(@(y) find(strcmp(y,processedOutput{2,paramCell}(:,1))>0),splitstr,'UniformOutput',false);
    ind= ~cellfun('isempty',oo);
    elementPos = find(ind >0);
    sel = cell2mat(oo(ind));
    splitstr(elementPos) = processedOutput{2,paramCell}(sel,3);
    
    mergeStr = cell(1,size(splitstr,2) + size(matchstr{1},2));
    mergeStr(1:2:end) = splitstr;
    mergeStr(2:2:end) = matchstr{1};
         
    processedOutput{2,fluxesCell}{k,2} = strjoin(mergeStr,'');
    
end

for k = 1:size(processedOutput{2,odeEqcell}(:,1))
    currentOdeEq = processedOutput{2,odeEqcell}(k,2);
    [matchstr,splitstr] = regexp(currentOdeEq,'\*|\+|(?<=\w|])\-(?=\w|[)|\/','match','split');
    % replacing parameters 
    oo = cellfun(@(y) find(strcmp(y,processedOutput{2,paramCell}(:,1))>0),splitstr{1},'UniformOutput',false);
    ind= ~cellfun('isempty',oo);
    elementPos = find(ind >0);
    sel = cell2mat(oo(ind));
    splitstr{1}(elementPos) = processedOutput{2,paramCell}(sel,3);
    
    mergeStr = cell(1,size(splitstr{1},2) + size(matchstr{1},2));
    mergeStr(1:2:end) = splitstr{1};
    mergeStr(2:2:end) = matchstr{1};
    
    processedOutput{2,odeEqcell}{k,2} = strjoin(mergeStr,'');
    
end

% print events
if sum(eventsCell) > 0
    fprintf(fid,'\n%%%% events %%%%%%\n');
    cellfun(@(x,y) fprintf(fid,'if (%s)\n %s = %s; \n end \n \n',x,y{1},y{2}),processedOutput{2,eventsCell}(:,1),processedOutput{2,eventsCell}(:,2),'UniformOutput',false);
end

% print fluxes
fprintf(fid,'\n%%%% fluxes %%%%%%\n');
cellfun(@(x,y) fprintf(fid,'%s=%s;\n',x,y),processedOutput{2,fluxesCell}(:,1),processedOutput{2,fluxesCell}(:,2),'UniformOutput',false);
fprintf(fid,'\n%%%% state equations %%%%%%\n');
cellfun(@(x,y) fprintf(fid,'d%s=%s;\n',x,y),processedOutput{2,odeEqcell}(:,3),processedOutput{2,odeEqcell}(:,2),'UniformOutput',false);
fprintf(fid,'dx = dx''; % return vector must be a column vector');
fclose(fid);
toc

% building the output structure 
dataOut.parameters = cellfun(@(x) str2num(x),processedOutput{2,3}(:,2));
% check each species has a corresponding ODE equation, one can add species
% to the model without reactions attach to it.
num = 1;
for k = 1:size(processedOutput{2,4}(:,1))
    ind = findStringInAList(processedOutput{2,odeEqcell}(:,1),processedOutput{2,4}(k,1));
    if ind ==0
        excludeInd(num) = k;
        num = num+1;
    end
end
selectedInitConds = ~ismember(1:size(processedOutput{2,4}(:,1)),excludeInd);
dataOut.initialValues = cellfun(@(x) str2num(x),processedOutput{2,4}(selectedInitConds,2));
dataOut.modelFile = [fileName '.m'];
dataOut.modelFilePath = filePath;
dataOut.modelFcn = fileName;
dataOut.speciesNames = processedOutput{2,odeEqcell}(:,1);
dataOut.parameterNames = processedOutput{2,paramCell}(:,1);
dataOut.reactionFluxes = processedOutput{4};



end

% parse expression type (A = B) lines 
function outCellArray = parseExpression(str,dividerStr)
    outCellArray = regexp(str,dividerStr,'split','once');
end
% parse Table like expressions 
function outCellArray = parseTable(str,dividerStr)
    expLoc = find(~cellfun('isempty', strfind(str,dividerStr)) >0);
    rr = regexp(str,'\t','split');
    outCellArray = rr(expLoc);
    
end


