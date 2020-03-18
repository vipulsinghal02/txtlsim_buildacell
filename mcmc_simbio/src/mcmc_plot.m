function mcmc_plot(marray, legendz, varargin)
% Plot the cornerplot of the parameter distributions, the autocorrelation
% function, and (optionally) the mcmc walker chains.
%
%
% INPUTS:
% marray: A numerical array of dumensions nParam x nWalkers x nSamp
% or nPoints x nParam. The walker chains are only plotted if the size is
% nParam x nWalkers x nSamp.
%
% legendz: A cell array of strings of parameter names, of size 1 x nParam.
%
% OPTIONAL NAME VALUE PAIR ARGUMENTS
%
%
% names: A cell of strings with labels for each parameter in m.
% ks: enable kernel smoothing density instead of histograms [false]
% support: a 2 x nParam matrix with upper and lower limits.
% scatter: show scatter plot instead of 2d-kernel density estimate
%           [true if #points<2000].
% fullmatrix: display upper corner of plotmatrix. [false]
% color: A color-theme for the plot. [.5 .5 .5].
% grid: show grid. [false].
%

% Copyright (c) 2018, Vipul Singhal, Caltech
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%
p = inputParser;
p.addParameter('nWalkers', [], @isnumeric);
% if empty then the length of the second dimension of marray is used.

p.addParameter('ks', false, @islogical);
p.addParameter('scatter', true, @islogical);
% p.addParameter('range', [99.5]);
p.addParameter('names', []);
p.addParameter('support', []); % if isempty, do nothing.
p.addParameter('ess', 100, @isnumeric);
p.addParameter('fullmatrix', false, @islogical);
p.addParameter('color', [.6 .35 .3]);
p.addParameter('grid', false, @islogical);
p.addParameter('fontsize', 18, @isnumeric);
p.addParameter('nbins', [], @isnumeric);
p.addParameter('transparency', 0.05, @isnumeric);
p.addParameter('savematlabfig', false, @islogical);
p.addParameter('savejpeg', false, @islogical);
p.addParameter('projdir', [], @ischar);
p.addParameter('tstamp', [], @ischar);
p.addParameter('plotChains', true, @islogical);

p.addParameter('extrafignamestring', [], @ischar);
p.parse(varargin{:})
p = p.Results;

% figure out the dimensions of the data array
if ndims(marray) == 3
    threeD = true;
elseif ndims(marray) == 2
    threeD = false;
else
    error('marray does not have the right number of dimensions.')
end

% Plot the corner plot

ecornerplot_vse(marray,'scatter', p.scatter,...
    'ks', p.ks,...
    'names', legendz,...
    'grid', p.grid,...
    'color', p.color,...
    'fullmatrix', p.fullmatrix,...
    'ks', p.ks,...
    'support', p.support, ...
    'transparency', p.transparency,...
    'fontsize', p.fontsize, ...
    'ess', p.ess);
if p.savematlabfig
    if isempty(p.projdir) || isempty(p.tstamp)
        warning('timestamp and project directory not specified. Nothing will be saved.')
    else
        
        specificproj = [p.projdir '/simdata_' p.tstamp];
        saveas(gcf, [specificproj '/cornerplot' p.tstamp p.extrafignamestring]);
        
    end
end

if p.savejpeg
    if isempty(p.projdir) || isempty(p.tstamp)
        warning('timestamp and project directory not specified. Nothing will be saved.')
    else
        specificproj = [p.projdir '/simdata_' p.tstamp];
        print(gcf, '-djpeg', '-r200', [specificproj '/cornerplot' p.tstamp p.extrafignamestring])
    end
end



% Plot the Ensemble Average Autocorrelation Function
%     figure
%     [C,lags,ESS]=eacorr(marray);

%     plot(lags,C,'.-',...
%         lags([1 end]),[0 0],'k');
%
%     grid on
%     xlabel('lags')
%     ylabel('autocorrelation');
%
%     text(lags(end),0,...
%         sprintf('Effective Sample Size (ESS): %.0f_ ',ceil(mean(ESS))),...
%         'verticalalignment','bottom',...
%         'horizontalalignment','right')


% Optionally plot the chains (if ndims(marray) == 3)
if ndims(marray)==3
    if p.plotChains
        if isempty(p.nWalkers)
            nW = size(marray, 2);
        else
            nW = p.nWalkers;
        end
        plotChains(marray, nW, legendz)
        
        
        if p.savematlabfig
            if isempty(p.projdir) || isempty(p.tstamp)
                warning('timestamp and project directory not specified. Nothing will be saved.')
            else
                
                specificproj = [p.projdir '/simdata_' p.tstamp];
                saveas(gcf, [specificproj '/trace' p.tstamp p.extrafignamestring]);
                
            end
        end
        
        if p.savejpeg
            if isempty(p.projdir) || isempty(p.tstamp)
                warning('timestamp and project directory not specified. Nothing will be saved.')
            else
                specificproj = [p.projdir '/simdata_' p.tstamp];
                print(gcf, '-djpeg', '-r200', [specificproj '/trace' p.tstamp p.extrafignamestring])
            end
        end
        
        
        
    end
    
    
    
end

end
