% protsynt81bis.m - MATLAB function for TXTL gene expression
% 
% Copyright (c) 2012 by California Institute of Technology
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%
% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
% 
% 3. Neither the name of the California Institute of Technology nor
%    the names of its contributors may be used to endorse or promote
%    products derived from this software without specific prior
%    written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CALTECH
% OR THE CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
% USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
% OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
% SUCH DAMAGE.
% 
% Original Author: Clare C. Chen
% Date: 5 September 2012

function dy=protsynt81bis(t,y, chgpar)
dy=zeros(2,1);

% model for the cell-free expression of 1 gene in our cell-free system
kcm = 0.05;         % rate of mRNA synthesis [s-1]
Km = 100;           % Michaelis constant for mRNA synthesis [nM]
S70t =5;            % total concentration of sigma70 [nM]
kdm=0.001;          % rate of mRNA degradation [s-1], 1/780 s-1
kcp=0.05;           % rate of protein synthesis [s-1]
Kp=500;             % Michaelis constant for protein synthesis [nM]
kdp=0;              % rate of protein degradation [nM s-1]
Et = 200;           % total concentration of active core E. coli RNAP [nM]
k70=0.26;           % dissociation constant sigma 70 and E. coli core RNAP [nM]
Lm=800;             % length of messenger RNA [bp]
Ld=868;             % length of DNA template that can be degraded [bp]
Cd=350;             % rate of DNA degradation by exonuclease [bp.s-1] 
Cm=6;               % rate of transcription [bp.s-1]
Cp=2;               % rate of translation [aa.s-1]
Rt = 2500;          % concentration of ribosomes [nM]
kf=0.002;           % rate of protein maturation [s-1], 1/480 s-1 for deGFP
kdx = 0.0000005;    % rate of degradation from enzymes other than RecBCD
Exot = 0.3;         % total exonuclease concentration [nM]
Kexo = 0.007;       % Michaelis constant for exonuclease degradation of DNA [nM] 
Gams = chgpar(1)*1000; % GamS supplied in micromolars, convert to nM
kgs = 0.71;         % dissociation constant for the binding of exonuclease to GamS[nM] 
ExoGams = Exot*Gams/(kgs + Gams); % concentration of exonuclease bound to GamS [nM]
Exotnew = Exot - ExoGams; % concentration of exonuclease not bound to GamS [nM]

%--------------------------------------------------------------------------
% Find the concentration of free RNA polymerase for a given DNA template
% concentration
%---------------------------------------------------------------------------
options = optimset('TolX', 1e-4);
e0=@(E0)E0+E0*S70t/(k70+E0)+E0*S70t*y(1)/(E0*S70t+Km*(k70+E0))*(1+kcm*Lm/Cm)-Et;
if sign(e0(0)*e0(Et)) > 0
    E0 = Et;
else
    E0=fzero(e0,[0,Et+1], options);
   
end

%---------------------------------------------------------------------------
% Find the free ribosome concentration for a given DNA template concentration
%---------------------------------------------------------------------------
r0=@(R0)R0+R0*y(2)/(R0+Kp)*(1+kcp*Lm/Cp)-Rt;
if sign(r0(0)*r0(Rt)) > 0
    R0 = Rt;
else
    R0=fzero(r0,[0,Rt + 1], options);
end

%---------------------------------------------------------------------------
% Find the concentration of P70 that is not bound to exonuclease at each time
% step for given DNA template concentrations.
%--------------------------------------------------------------------------
p70f = @(P70f)(P70f/Kexo)*(Exotnew - (y(1) - P70f)) - (y(1) - P70f);
if sign(p70f(0)*p70f(y(1))) > 0
    P70f = y(1);
else
    P70f = fzero(p70f, [0, y(1)]);
end

%--------------------------------------------------------------------------
% Find the concentration of P70 that is bound to exonuclease at each time
% step for each of the DNA template concentrations.
%-------------------------------------------------------------------------

P70b = y(1) - P70f;

%-------------------------------------------------------------------------
% define ODEs for system
%-------------------------------------------------------------------------
dy(1) = -(Cd/Ld)*P70b-kdx*y(1);			% DNA template concentration
dy(2) = kcm*y(1)*E0*S70t/(E0*S70t+Km*(E0+k70))-kdm*y(2);  % mRNA concentration
dy(3) = kcp*y(2)*R0/(Kp+R0)-kdp-kf*y(3);	% non mature GFP concentration
dy(4) = kf*y(3)-kdp;				% mature GFP concentration
