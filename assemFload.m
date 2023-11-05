%**************************************************************************
% Assemble the pre-load depending on Ut in equation of Utt
%**************************************************************************
% DESCRIPTION
% Assemble the pre-load Fk depending on ut_ij in equation of utt_ijk, k=1,2
% Expression (with eta test function) : 
% """ int_Omega ( gamma*C*( ut x_s e_k ).(strain( eta )) """
%   
% CALLED in :
% ... solve2.m
%
% INPUT
%  - mesh:     pdetool mesh struct
%  - p,t:      pdetool p,t
%  - tgamma2:  contrast in center of triangles
%  - la0, mu0: pdetool coeficiets struct
%  - Ut:       first localisator ut_ij
%
% OUTPUT
%  - [F1,F2]: pre-loads F1 (k=1) and F2 (k=2)
%
% HISTORY
% A.A. Novotny, V. Calisti : 03/2020 : code implementation.
%**************************************************************************

function [F1,F2] = assemFload(mesh,p,t,tgamma2,la0,mu0,Ut)

tgamma2 = tgamma2.*[mesh.area;mesh.area];
[dsx,dsy] = sfderivatives(mesh);
dsx = dsx';
dsy = dsy';
% Utg(1,:) = utx ; Utg(2,:) = uty 
Utg = pdeintrp(p,t,Ut); 
% multiply by gamma*area here for later integration
Utg = tgamma2.*Utg;
% usefull coefficient
ct0  = 2*mu0 + la0;

% for j=1:np
% for F1 we calculate : F1_(j) and F1_(j+np) 
%   =ct0*Utg(1,:)*Dx(phix) + la0*Utg(1,:)*Dy(phiy)           ...
%                           + mu0*Utg(2,:)*(Dx(phiy)+Dy(phix))  ;
% for F2 we calculate : F2_(j) and F2_(j+np)
%   =la0*Utg(2,:)*Dx(phix) + ct0*Utg(2,:)*Dy(phiy)           ...
%                           + mu0*Utg(1,:)*(Dx(phiy)+Dy(phix))  ;
% where we have :   phi = [ phix phiy ]  (len=2*np)
%     * for j     : phi = [ PHIj   0   ]  ---> Fx
%     * for j+np  : phi = [  0   PHIj  ]  ---> Fy
% with PHIj = [0 ... 0 1 0 ... 0]' jeme coord, (len=np) ; Dx=dsx ; Dy=dsy

%--------F1
% F1x
F1x = ct0*dsx*(Utg(1,:)') + mu0*dsy*(Utg(2,:)')  ;
% F1y
F1y = la0*dsy*(Utg(1,:)') + mu0*dsx*(Utg(2,:)')  ;
% F1
F1 = [F1x ; F1y];
%--------F2
% F2x
F2x = la0*dsx*(Utg(2,:)') + mu0*dsy*(Utg(1,:)')  ;
% F2y
F2y = ct0*dsy*(Utg(2,:)') + mu0*dsx*(Utg(1,:)')  ;
% F2
F2 = [F2x ; F2y];

end