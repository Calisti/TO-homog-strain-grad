%***********************************************************************
% Finite element stress calculation
%***********************************************************************
% DESCRIPTION 
% Voigt convention :
% s = 2mu*e + la*tr(e)I  where e is the strain tensor
% s = [ s11 ;
%       s22 ;
%       s12  ]
%***********************************************************************

function s=getstress(u,p,t,la,mu)
e  = getstrain(u,p,t); 
id = [1 1 0]';
s  = la*id*(e(1,:) + e(2,:)) + 2*mu*e;
end