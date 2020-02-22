% Zoltan A Tuza July 2013
%
% This function handles biotek plate reader data in csv format. The delimiter in
% the file must be given. The third optional argument is the well of the negative
% control. By default the 2nd well is treated as negative control.
%
%
% Copyright (c) 2013 by California Institute of Technology
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

function experimentData = wholeExpfileReader(fileName,delimiter,varargin)

data = csvimport(fileName,'delimiter',delimiter);


if nargin > 2
    noBGwell = varargin{1};
else
    noBGwell = 2;
end

experimentData.Bgwell = noBGwell;
experimentData.FileName = fileName;


%%

%find the string with runtime and Interval info  
res = cellfun(@(x) sscanf(x,'Runtime %d:%d:%d %*s Interval %d:%d:%d, %d Reads '),data(1:100,:),'UniformOutput',false);

res2 = cellfun(@isempty,res);
d = cell2mat(res(find(res2 == 0)));

experimentData.runTime = [3600 60 1]*d(1:3); % stored in seconds
experimentData.readInterval = [3600 60 1]*d(4:6); % stored in seconds
experimentData.numOfReads = d(7);
experimentData.numOfValidReads = d(7);



res = cellfun(@(x) strcmp(x,'Read'),data(1:100,:),'UniformOutput',false);
ind = find(cell2mat(res) == 1);


for k = 1:size(ind,1)
    experimentData.channels(k,1) = data(ind(k),2);
    experimentData.channels(k,1) = data(ind(k),2);
    exemVec = sscanf(data{ind(k)+4,2},'Excitation: %d, Emission: %d');
    experimentData.channels(k,2) = {exemVec'};
    gain = sscanf(data{ind(k)+5,2},'Optics: %*s Gain: %d');
    experimentData.channels(k,3) = {gain};
    
    % mine out data matrices
    res = cellfun(@(x) strcmp(x,sprintf('%s:%d,%d',experimentData.channels{k,1},experimentData.channels{k,2}) ),data,'UniformOutput',false);
    dataHeaderInd = find(cell2mat(res) == 1);
    % determine the number of data columns
    wellNum = sum(~cellfun(@isempty,data(dataHeaderInd+2,4:end)));
    experimentData.wellNames = data(dataHeaderInd+2,(1:wellNum)+3);
    % time vector
    tVecStr = data(dataHeaderInd+3:(dataHeaderInd+3+experimentData.numOfReads)-1,2);
    experimentData.t_vec = cellfun(@(x) [3600 60 1]*sscanf(x,'%d:%d:%d'),tVecStr);
    % handle canceled runs
    tmpData = data(dataHeaderInd+3:(dataHeaderInd+3+experimentData.numOfReads)-1,(1:wellNum)+3);
    tmpData = cellfun(@(x) str2num(x),tmpData,'UniformOutput',false);
    emptyIndex = cellfun(@isempty,tmpData); % Find indices of empty cells
    experimentData.numOfValidReads = experimentData.numOfReads - sum(emptyIndex(:,1));
    % put zeros to empity cells -> needed for cell2mat conversion
    tmpData(emptyIndex) = {0};
    experimentData.Data(:,:,k) = cell2mat(tmpData);
    
    experimentData.noBg(:,:,k) = experimentData.Data(:,:,k) - repmat(experimentData.Data(:,noBGwell,k),1,wellNum);
    
end

% experiment can be canceled between channel reads, in that case the
% maximal number of full cycle is treated as numOfValidReads
experimentData.noBg = experimentData.noBg(1:experimentData.numOfValidReads,:,:);
experimentData.t_vec = experimentData.t_vec(1:experimentData.numOfValidReads);

% removing MG channel background signal by curve fitting
% search for mg channel
channel_index = find(cellfun(@(x) strcmpi(x,'MG'),experimentData.channels(:,1)) > 0);
if ~isempty(channel_index)
    % on the biotek plate reader there are two distinct processes affecting
    % background -> two gaussian function covers that.
    %f = fittype('gauss2');
    % new approach: smoothing spline
    f2 = fittype( 'smoothingspline' );
    opts = fitoptions( f2 );
    opts.SmoothingParam = 1.353e-6;
    [fitinfo,gof] = fit(experimentData.t_vec/60,experimentData.Data(1:experimentData.numOfValidReads,experimentData.Bgwell,channel_index),f2,opts);
    bgAdjusted = experimentData.Data(1:experimentData.numOfValidReads,:,channel_index) - repmat(fitinfo(experimentData.t_vec/60),1,size(experimentData.Data,2));
    experimentData.noBg(:,:,channel_index) = bgAdjusted;
    experimentData.MGBgFitinfo = fitinfo;
    experimentData.MGBgGof = gof;
end

% calculate endTime, expression rate
for k=1:size(experimentData.channels,1)
    dataColumnNum = size(experimentData.noBg(:,:,k),2);
    rate_x = num2cell(repmat(experimentData.t_vec,1,dataColumnNum),[1 dataColumnNum]);
    rate_y = num2cell(experimentData.noBg(:,:,k),[1 dataColumnNum]);

    [experimentData.rate(:,k) experimentData.rateFitInfo(:,k)] = cellfun(@(x,y) expression_rate(x,y),rate_x,rate_y,'UniformOutput',false);
    
    [experimentData.endTimeIndex(:,k) experimentData.endTime(:,k)] = cellfun(@(x,y) expression_endtime(x,y),rate_x,rate_y,'UniformOutput',false);
    
end

% convert back matrices
experimentData.rate = cell2mat(experimentData.rate);
experimentData.endTimeIndex = cell2mat(experimentData.endTimeIndex);
experimentData.endTime = cell2mat(experimentData.endTime);





