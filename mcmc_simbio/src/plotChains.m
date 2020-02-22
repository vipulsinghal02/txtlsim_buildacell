function plotChains(m, nW, legends, varargin)
% Plot markov chain trajectories. 
% m is a nParam x nWalkers x numSamples
% nW is the number of walkers to plot. 
% legs is legends
p = inputParser;
p.addParameter('Visible', 'on', @ischar)
p.parse(varargin{:});
p=p.Results;

[nParam, nWalkers, nSamples] = size(m);
[n1 n2] = twofactors(nParam);

if isprime(nParam)
    n1 = ceil(nParam/5); % 5 columns
    n2 = 5;
end

wix = unique(ceil(rand(nW, 1)*nWalkers));
m = m(:,wix, :);
figure('Visible', p.Visible)
ss = get(0, 'screensize');
set(gcf, 'Position', [ss(3)*(1-1/1.1) ss(4)*(1-1/1.15) ss(3)/1.1 ss(4)/1.15]);
for i = 1:nParam
    subplot(n1, n2, i)

    for j = 1:length(wix)
        plot(1:nSamples, squeeze(m(i, j, :)),...
            'LineWidth', 0.1,...
            'color', [0.2 0.7 0.1].^2)
        hold on       
    end
    title(legends{i}, 'fontsize', 18);
    xlabel('Iteration', 'fontsize', 16);
    ylabel('log-value', 'fontsize', 16);
    ax = gca;
    ax.FontSize = 16;
    
    
end

end

