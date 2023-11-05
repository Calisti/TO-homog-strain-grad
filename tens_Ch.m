%**************************************************************************
% Computes the homogenised tensor Ch
%**************************************************************************
% DESCRIPTION
% Given the canonical stress field, computes the 
% homogenised elasticity tensor Chom
% Chom is symetric
%            kl=11    kl=22    kl=12
% Chom = [ Ch_1111  Ch_1122  Ch_1112  ;   ij = 11
%             *     Ch_2222  Ch_2212  ;   ij = 22
%             *         *    Ch_1212   ]  ij = 12
%
%   ch11 ch12 ch13
%   ch21 ch22 ch33
%   ch31 ch32 ch33
%
% INPUT
%  - mesh:        pdetool mesh struct
%  - U:           micro displacement
%  - matprop:     material properties struct
%  - psi:         level-set function
%  - cell_volume: volume of the REV
%
% OUTPUT
%  - Ch:   homogenized elasticity tensor Chom 
%
% HISTORY
% A.A. Novotny, V. Calisti      03/2020: code updating.
%***********************************************************************

function Ch = tens_Ch(mesh,U,matprop,psi,cell_volume)

    gamma = matprop.gamma; la0 = matprop.la0; mu0 = matprop.mu0;
    % element characteristic function: 1 if psi > 0 and gama if psi <= 0
    tgamma = pdeintrp(mesh.p,mesh.t,(psi<0) + gamma*(psi>=0));
    aux = tgamma.*mesh.area; tgamma3 = [aux;aux;aux];
    ch = zeros(3,3);
    for i=1:3
        % i = ij
        u=U(:,i);
        %nominal stress
        s=getstress(u,mesh.p,mesh.t,la0,mu0); 
        aux=tgamma3.*s; ch(i,:)=sum(aux,2);
    end
    aux = ch/cell_volume;
    
    % e = (e11, e22, e12)
    ch=[[aux(1,1),aux(1,2),aux(1,3)];
        [aux(2,1),aux(2,2),aux(2,3)];
        [aux(3,1),aux(3,2),aux(3,3)]];
    
    Ch.tens = ch;
%     [Ch.eigV,Ch.eigv]  = eig(ch);
    Ch.nm   = norm(ch);
    
end
