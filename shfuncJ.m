%***********************************************************************
% Computation of the value for each shape functional
%***********************************************************************
% DESCRIPTION
% Computes the value for each shape functional
% 
% INPUT
%  - cost:  option for the choice of the shape functional
%  - Ch,Dh: homogenised elasticity tensor
%
% OUTPUT
%  - J: value of the choosen shape functional
%
% HISTORY
% S. Amstutz & A.A. Novotny      06/2009: code implementation.
% S.M. Giusti                    03/2011: code updating. 
% V. Calisti                     2021: code updating. 
%***********************************************************************

function J = shfuncJ(cost,Ch,Dh)

Ch = Ch.tens ;
Dh = Dh.tens ;
    
FunJ = cost.FunJ;
Idx  = cost.Idx;
IC   = Idx.IC   ; IS   = Idx.IS   ; IDm = Idx.IDm ;
Sh = 0;
if ~ isempty(IS)
    Sh = inv(Ch);
end
J = FunJ([ Ch(IC) ; Sh(IS) ; Dh(IDm) ]);

end
