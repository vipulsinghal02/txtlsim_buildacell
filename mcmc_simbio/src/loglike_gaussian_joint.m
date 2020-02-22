function [llh] = loglike_gaussian_joint(jointlogp, stdev, ICarray,measuredspidx, tspan, ...
	dataArray1,dataArray2, genemodel, lognormvec)
% 
% loglike_gaussian Gaussian log likelihood
%   Compute the log likelihood of the data given the model parameters and
%   the standard deviation. 
%
% Data is input as a array that is  N x M x nIC. N is number of time points in tspan, 
% M is the number of species species
% 
% one param vector, used two populate two runs. two data sets. 
logp1 = jointlogp([1:4 7:8]);

mV = mean(dataArray1(:,measuredspidx,:),1);
meanVals = mean(mV,3); %meanVals is a row vec.

wt = sum(meanVals)./meanVals; %hight mean = lower wt
relWt = wt/sum(wt);
CONC_temp = zeros(size(dataArray1));

nIC = size(ICarray,1); % number of different initial conditions. 
% initial conditions are row vectors. Rows are different sets of initial conditions.

for i = 1:nIC
[~,CONC_temp(:,:,i)] = genemodel(logp1, ICarray(i,:), tspan);
end

residuals = dataArray1(:,measuredspidx,:) - CONC_temp(:,measuredspidx,:);
resi = bsxfun(@times, residuals, relWt);

res1 = resi(:);

logp2 = jointlogp([1:2 5:8]);
mV = mean(dataArray2(:,measuredspidx,:),1);
meanVals = mean(mV,3); %meanVals is a row vec.

wt = sum(meanVals)./meanVals; %hight mean = lower wt
relWt = wt/sum(wt);
CONC_temp = zeros(size(dataArray2));


nIC = size(ICarray,1); % number of different initial conditions. 
% initial conditions are row vectors. Rows are different sets of initial conditions.

for i = 1:nIC
[~,CONC_temp(:,:,i)] = genemodel(logp2, ICarray(i,:), tspan);
end

residuals = dataArray2(:,measuredspidx,:) - CONC_temp(:,measuredspidx,:);
resi = bsxfun(@times, residuals, relWt);

res2 = resi(:);

res = [res1; res2];
llh = sum(lognormvec(res, stdev));



end

