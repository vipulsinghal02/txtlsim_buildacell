function marray = mcmc_get_walkers(tstamp, simID, projdir)
% get mcmc walker samples from saved data.
%
% - tstamp: a cell array of strings (N x 1). Each string
% corresponds to a time stamp of the simulation where the data is stored.
% Single strings should still be encapsulated in cells.
%
% Example:
%
% {'20180117_173057', '20180117_164124', '20180117_173057'}
%
% OR
%
% {'20180117_173057'} <--- singe timestamp still incapsulated in the cell.
%
%
% - simID. A cell array of size N x 1. The ith cell of this array can contain
% either a vector of numerical indices (of dimension Mi x 1),
% or a cell array of dimension Mi x 1, where each cell has a string specifying
% the simulation associated with the ith
% timestamp in the array tstamp. Single strings, or a single integer vector array
% should still be encapsulated in a cell.
%
%
%
% - projdir is the directory where the folders simdata_tstamp{i} are stored.

% ------------------------------------------

% Copyright (c) 2018, Vipul Singhal, Caltech
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

if iscell(tstamp)
    if all(cellfun(@iscell, simID))
        mcat = retrievemcat(projdir, tstamp, simID);
    elseif all(cellfun(@isnumeric, simID))
        % convert simID to a cell array of cell arrays of strings. (from a
        % cell array of numbrical vector arrays)
        convertedsID = cellfun(...
            @(numarray) arrayfun(@num2str,numarray, 'UniformOutput', false),...
            simID,...
            'UniformOutput', false);
        mcat = retrievemcat(projdir, tstamp, convertedsID);
    else
        % simID not a cell array when the timestamps are. Needs to be a cell array.
        error('simID must be a cell array of cell arrays of strings, or a cell array of integer vector arrays.')
    end
else
    error('tstamp must be a cell array of strings.');
end
marray = mcat;
clear mcat;
clear m;

end

% Define local Functions.

function mcat = retrievemcat(projdir, tstamp, simID)
% simdata must be a cell of strings of simulation IDs (numbers that identify the simulation number for a p
% particular timestamp. )

load([projdir '/simdata_' tstamp{1} '/mcmc' tstamp{1} '_ID' simID{1}{1} '.mat'], 'm');
mcat = m;
for j = 2:length(simID{1})
    % load , check number of parameters and walkers, and concatenate.
    load([projdir '/simdata_' tstamp{1} '/mcmc' tstamp{1} '_ID' simID{1}{j} '.mat'], 'm');
    % check dimensions
    mcat = cat(3, mcat, m);
end
for i = 2:length(tstamp)
    load([projdir '/simdata_' tstamp{i} '/mcmc' tstamp{i} '_ID' simID{i}{1} '.mat'], 'm');
    if ~isequal(size(mcat, 2), size(m, 2)) % # of walkers are different!
        warning('The number of walkers in arrays to be concatenated do not match; rearranging to match. ')
        rearrangementFlag = true;
        mtemp = m(:, :);
        for j = 2:length(simID{i})
            load([projdir '/simdata_' tstamp{i} '/mcmc' tstamp{i} '_ID' simID{i}{j} '.mat'], 'm');
            mtemp = [mtemp m(:,:)];
        end
        nW = size(mcat, 2);
        ntemp = size(mtemp, 2);
        np = size(mcat, 1);
        nplanes = floor(ntemp/nW);
        m2temp = mtemp(:,1:nplanes*nW);
        m3 = reshape(m2temp, np, nW, nplanes);
        mcat = cat(3, mcat, m3);
    else
        mcat = cat(3, mcat, m);
    end
    % massive bugfix (feb5, 2019): forgot to add the code for the
    % second ID onwards.
    
    for j = 2:length(simID{i})
        % load , check number of parameters and walkers, and concatenate.
        load([projdir '/simdata_' tstamp{i} '/mcmc' tstamp{i} '_ID' simID{i}{j} '.mat'], 'm');
        if ~isequal(size(mcat, 2), size(m, 2)) % # of walkers are different!
            warning('The number of walkers in arrays to be concatenated do not match; rearranging to match. ')
            rearrangementFlag = true;
            mtemp = m(:, :);
            for kk = 2:length(simID{i})
                load([projdir '/simdata_' tstamp{i} '/mcmc' tstamp{i} '_ID' simID{i}{kk} '.mat'], 'm');
                mtemp = [mtemp m(:,:)];
            end
            nW = size(mcat, 2);
            ntemp = size(mtemp, 2);
            np = size(mcat, 1);
            nplanes = floor(ntemp/nW);
            m2temp = mtemp(:,1:nplanes*nW);
            m3 = reshape(m2temp, np, nW, nplanes);
            mcat = cat(3, mcat, m3);
        else
            mcat = cat(3, mcat, m);
        end
    end
    
end
end

