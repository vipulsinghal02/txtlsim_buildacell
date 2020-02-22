%% data importer
% Zoltan A Tuza Feb 2013
%
% This function handles plate reader data in csv format. The delimiter in
% the file has to be ';'. Please, save your file accordingly. 
%
% Victor plate reader data should be ordered by valves before import
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

function experim =  data_importer(filename,plate_reader,varargin)

if ~exist(filename,'file')
    error('%s not exist!',filename)
end

experim = importdata(filename,';',1);


% check in the biotek case temparature column is included, or not
temp = strfind(experim.textdata(1,2),'T°');
   if ~isempty(temp{1})
    noBgValve = 3;
    startValve = 2;
   else 
    noBgValve = 2;
    startValve = 1;
   end 


% If not define the 2nd valve is treated as negative control valve.
if ~isempty(varargin)
   noBgValve =  varargin{1};
   startValve = 1;
end

%%%%%%%%%%%%%%%%%%% Import Victor data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(plate_reader,'victor')
    repeats = str2num(experim.textdata{end,2});
    numOfValves = size(experim.data,1)/repeats;
    experim.t_vec = cellfun(@(x) [3600 60 1 0]*sscanf(x,'%d:%d:%d.%d'),experim.textdata(2:repeats+1,5));
    
    experim.valves = reshape(experim.data,repeats,numOfValves);
    
    experim.noBg = experim.valves - repmat(experim.valves(:,noBgValve),1,numOfValves);
 
%%%%%%%%%%%%%%%%%%% biotek data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif strcmp(plate_reader,'biotek')
    experim.t_vec = cellfun(@(x) [3600 60 1]*sscanf(x,'%d:%d:%d'),experim.textdata(2:end,1));
    
    experim.noBg = experim.data(:,startValve:end) - repmat(experim.data(:,noBgValve),1,size(experim.data(:,startValve:end),2));
%%%%%%%%%%%%%%%%%%% specs data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif strcmp(plate_reader,'specs')
    
else
    error('currently only ''victor '' and ''biotek'' options are supported ')
    
end


end





