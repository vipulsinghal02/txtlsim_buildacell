% All parameter values were taken from 
% Nature. 2000 Jan 20;403(6767):335-8.

% This model was developed by members of the 2003 Synthetic Biology Class
% on Engineered Blinkers

% set the max transcription rate in transcripts per second
k_transcription_cI = 0.5;
k_transcription_lacI = 0.5;
k_transcription_tetR = 0.5;

% set the leakage transcription rate (ie transcription rate if
% promoter region bound by repressor) in transcripts per second
k_transcription_leakage = 5e-4;

% Set the mRNA and protein degradation rates (per second)
mRNA_half_life = 120; % in seconds
k_mRNA_degradation = log(2)/mRNA_half_life;
protein_half_life = 600; % in seconds
k_protein_degradation = log(2)/protein_half_life;

% proteins per transcript lifespan
translation_efficiency = 20; 

average_mRNA_lifespan = 1/k_mRNA_degradation; % in seconds

% proteins per transcript per sec
k_translation = translation_efficiency/average_mRNA_lifespan; 

% set the Hill coefficients of the repressors
n_tetR = 2;
n_cI = 2;
n_lacI = 2;

% Set the dissociation constant for the repressors to their target promoters
% in per molecule per second
KM_tetR = 40;
KM_cI = 40;
KM_lacI = 40;

