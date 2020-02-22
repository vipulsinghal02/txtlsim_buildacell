function mcat = catMC(datafiles)
%catMC Take Markov Chains from GWMCMC and concat them

load(datafiles{1}, 'm');
mcat = m;
for i = 2:length(datafiles)
    load(datafiles{i}, 'm');
    mcat = cat(3, mcat, m);
end

