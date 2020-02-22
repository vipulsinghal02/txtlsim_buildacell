delete(gcp('nocreate'))
parpool(48)
primeNumbers = primes(uint64(2^23));
compositeNumbers = primeNumbers.*primeNumbers(randperm(numel(primeNumbers)));
factors = zeros(numel(primeNumbers),2);
numWorkers = [12 18 24 32 36 46 47 48];
tLocal = zeros(size(numWorkers));
for w = 1:numel(numWorkers)
    tic;
    parfor (idx = 1:numel(compositeNumbers), numWorkers(w))
        factors(idx,:) = factor(compositeNumbers(idx));
    end
    tLocal(w) = toc;
end
numWorkers2 = [numWorkers]
tLocal2 = [tLocal]
f = figure;
speedup = tLocal2(1)./tLocal2;
plot(numWorkers2, speedup);
title('Speedup with the number of workers');
xlabel('Number of workers');
xticks(numWorkers2);
ylabel('Speedup');
%%
f = figure;

plot(numWorkers2, tLocal2);
title('time taken vs the number of workers');
xlabel('Number of workers');
xticks(numWorkers2);
ylabel('Time Taken');