function [bds, mx, mp, BW] = MPinterval(vec, varargin)
%[bds, mx, mp ] = MPinterval(vec, varargin) compute the smallest interval containing maximum probability estimate for a given
%probability mass
% OUTPUTS
% bds are the lower and upper bounds of the interval.
% mx is the sample at which the maximum probability occurs
% mp is the probability
% INPUTS
% vec is the vector of samples
% mass is a number from 0 to 1. 0 returns just the map estimate, and bds =
% [mx mx], 1 returns the entire interval
% nBins is the number of bins to use.

p = inputParser;
p.addParameter('mass',0.2,@isnumeric);
npoints = round(length(vec)/1000);
p.addParameter('npoints',npoints,@isnumeric);
p.addParameter('bw',[],@isnumeric);
p.parse(varargin{:});
p=p.Results;

% compute density, mp, and interval

if isempty(p.bw)
    [F,XI,BW]=ksdensity(vec, 'npoints', p.npoints);
else
    [F,XI]=ksdensity(vec, 'npoints', p.npoints, 'bandwidth', p.bw);
    BW = p.bw;
end

[mp, mxi] = max(F); % mxi is the index of the bin at which the max occurs
mx = XI(mxi);
currmass = 0;
currL = mxi;
currR = mxi;
while (currmass < p.mass) && currL>1 && currR<length(F)
    if F(currL) >= F(currR)
        currL = currL - 1 ;
    else
        currR =currR + 1;
    end
    newIdx = currL:currR;
    try
        currmass = trapz(XI(newIdx), F(newIdx));
    catch
        disp('error...')
    end
    
end

bds = [XI(currL) XI(currR)];
% note that the limits of the XI need not coincide with the limits of the
% the samples, since XI is the support of the gaussian fitted to the
% samples.
end

