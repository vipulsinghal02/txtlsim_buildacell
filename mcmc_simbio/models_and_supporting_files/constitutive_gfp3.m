function d_sp = constitutive_gfp3( t, sp, logp)
% constitutive_gfp3: just dna and enzyme model of gfp expression
% D_G + P <-> D_G__P  (kfP, krP)
% D_G__P -> D_G + P + G 
% 
% no conservation laws. 
% ESP: k_c 
% CSP: k_fP, k_rP
% Env specific species (ESSP): P
% Dosed species: D_G
% initial condition 0 species: D_G__P, G
%
% Use thise for setting up teh variables for the function
% log_likelihood_sharedCSP
% espix = 1;
% cspix = 2:3;
% esspix = 2; 
% pmap_calib = {espix, cspix, esspix};
% nSp_calib = 4; 


p = exp(logp);

% ESP
kc = p(1);

% CSP
kfP = p(2);
krP = p(3);

% Species 
D_G = sp(1); % dosed
P = sp(2); % ESSP
D_G__P = sp(3);
G = sp(4); % measured species. idMS = 4 in the function log_likelihood_sharedCSP


% ODEs
dD_G = -kfP*D_G*P + (krP+kc)*D_G__P;
dP = -kfP*D_G*P + (krP+kc)*D_G__P;
dD_G__P =  +kfP*D_G*P - (krP+kc)*D_G__P;
dG = kc*D_G__P;


d_sp = [dD_G; dP; dD_G__P; dG];

end

