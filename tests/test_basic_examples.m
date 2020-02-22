%% test file for basic examples
% Zoltan A Tuza Aug 2013
%
clearvars
close all
clc
%% gene expression
disp('testing gene expression')
geneexpr
GFP_ind = findspecies(Mobj,'protein deGFP*');
if GFP_ind ~= 0
    disp('GFP found');
    if sum(x_ode(:,GFP_ind)) >0
        disp('gene expression: PASSED');
    else
        disp('gene expression: FAILED');
    end
else
    disp('gene expression: FAILED');
end
clearvars
close all

%% negautoreg
disp('testing negative autoregulation')
negautoreg
tetR_ind = findspecies(Mobj,'protein tetRdimer');
if tetR_ind ~= 0
    disp('tetR found');
    if sum(x_ode(:,tetR_ind)) >0
        disp('negautoreg: PASSED');
    else
        disp('negautoreg: FAILED');
    end
else
    disp('negautoreg: FAILED');
end
clearvars
close all

%% combinatorial promoter
disp('testing combinatorial promoter')
combinatorial_promoter
GFP_ind = findspecies(well_a1,'protein deGFP-lva-terminator*');
if GFP_ind ~= 0
    disp('tetR found');
    if sum(x_ode(:,GFP_ind)) >0
        disp('combinatorial promoter: PASSED');
    else
        disp('combinatorial promoter: FAILED');
    end
else
    disp('combinatorial promoter: FAILED');
end
clearvars
close all

