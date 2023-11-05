%***********************************************************************
% Finite element strain calculation
%***********************************************************************
% DESCRIPTION 
% Voigt convention :
% e = [        dux/dx         ;         e = [ e11 ;
%              duy/dy         ;               e22 ;
%      (1/2)(dux/dy + duy/dx)  ]              e12  ]
%***********************************************************************

function e=getstrain(u,p,t)
[ux,uy] = pdegrad(p,t,u);
e       = [ux(1,:) ; uy(2,:) ; (ux(2,:)+uy(1,:))/2];
end