% Simplified model of the repressilator

% This function implements the basic model of the repressilator
% All parameter values were taken from 
% Nature. 2000 Jan 20;403(6767):335-8.

% This model was developed by members of the 2003 Synthetic Biology Class
% on Engineered Blinkers

function dydt = repressilator(t, y)

%store the state variables under more meaningful names
mRNA_cI = y(1);
mRNA_lacI = y(2);
mRNA_tetR = y(3);
protein_cI = y(4);
protein_lacI = y(5);
protein_tetR = y(6);

% set the parameter values
genereg_params;

% the differential equations governing the state variables:
% mRNA concentration = transcription given repressor concentration - 
% mRNA degradation + transcription leakage
dydt(1,:) = k_transcription_cI/(1 + (protein_tetR / KM_tetR)^n_tetR) - ...
  k_mRNA_degradation*mRNA_cI + k_transcription_leakage;
dydt(2,:) = k_transcription_lacI/(1 + (protein_cI / KM_cI)^n_cI) - ...
  k_mRNA_degradation*mRNA_lacI + k_transcription_leakage;
dydt(3,:) = k_transcription_tetR/(1 + (protein_lacI / KM_lacI)^n_lacI) - ...
  k_mRNA_degradation*mRNA_tetR + k_transcription_leakage;

% protein concentration = translation - protein degradation
dydt(4,:) = k_translation*mRNA_cI - k_protein_degradation*protein_cI;
dydt(5,:) = k_translation*mRNA_lacI - k_protein_degradation*protein_lacI;
dydt(6,:) = k_translation*mRNA_tetR - k_protein_degradation*protein_tetR;
