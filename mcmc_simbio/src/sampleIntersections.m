function [Ikeep, bds, mx, mp] = sampleIntersections(m, colIdx, varargin)
%sampleIntersections Indices of the MCMC samples to use for plotting
%   generate sample intersections using either custom bounds, percentile
%   bounds, or MAP estimate mass bounds
% m is numSamples x numParameters matrix
% colIdx is the indices of the parameters to to use

p = inputParser;
p.addParameter('mode','MAPmass',@ischar); % other modes: percentile, custom
p.addParameter('npoints',100,@isnumeric);
p.addParameter('bw',[],@isnumeric);
p.addParameter('mass',0.2,@isnumeric);

p.parse(varargin{:});
p=p.Results;

if strcmp(p.mode, 'MAPmass')
    bds = zeros(2, length(colIdx));
    mx = zeros(1, length(colIdx));
    mp = zeros(1, length(colIdx));
    Ikeep = (1:size(m, 1))';
    for i = 1:length(colIdx)
        if isempty(p.bw)
            [bds(:,i), mx(i), mp(i) ] = MPinterval(m(:,colIdx(i)), 'npoints', p.npoints,...
                'mass', p.mass);
        else
            [bds(:,i), mx(i), mp(i) ] = MPinterval(m(:,colIdx(i)), 'npoints', p.npoints,...
                'mass', p.mass, 'bw', p.bw);
            BW(i) = p.bw;
        end
       

        Inew = find(m(:,colIdx(i))>bds(1,i) & m(:,colIdx(i))<bds(2,i));
        Ikeep = intersect(Ikeep, Inew);
    end
end





end

