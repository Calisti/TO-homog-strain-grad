%**************************************************************************
% Updates pdetool linear form for periodic bc
%**************************************************************************
% DESCRIPTION
% Reduce the second member tens_F of linear equation KX=tens_F for periodic bc
% Original problem: 
% [ Kii Ki+ Ki- Kic;  *  [ Xi;   =  F  = [ Fi;
%   K+i K++ K+- K+c;       X+;             F+;
%   K-i K-+ K-- K-c;       X-;             F-;
%   Kci Kc+ Kc- Kcc ]      Xc ]            Fc ]
% Reduced problem:
% [   Kii    ( Ki+ +Ki-);   * [Xi;   = Freduced = [ Fi     ;
%   ( K+i    ( K++ +K+-        X+ ]                 F+ + F- ]
%    +K-i )   +K-+ +K--) ]      
%
% CALLED in :
% ...solve2.m
%
% INPUT
%  - F:    pdetool load vector
%  - bc:   extra boundary conditions struct
%
% OUTPUT
%  - Fr: updated pdetool load vector for periodic bc
%
% HISTORY
% A.A. Novotny, V. Calisti   03/2020: code implementation.
%**************************************************************************

function Fr = reducedF(F,bc)

dofs = bc.dofs;
ni   = size(dofs.i, 1);
np   = size(dofs.p, 1);

Fi  = F(dofs.i(:));
Fp  = F(dofs.p(:));
Fm  = F(dofs.m(:));
Fr  = zeros(ni+np, 1);
Fr  = [ Fi ; Fp + Fm ];
clear Fi Fp Fm
    
end
