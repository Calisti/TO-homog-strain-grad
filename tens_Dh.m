%**************************************************************************
% Computes the homogenised tensor D^hom
%**************************************************************************
% DESCRIPTION
% The coefficients of the homogenised tensor D^hom are saved thi way
% dh = (D_ijk,pqr)  Eij,k Epq,r
% = D111.111 D111.112   D111.221 D111.222   D111.121 D111.122
%   D112.111 D112.112   D112.221 D112.222   D112.121 D112.122
%
%   D221.111 D221.112   D221.221 D221.222   D221.121 D221.122
%   D222.111 D222.112   D222.221 D222.222   D222.121 D222.122
%
%   D121.111 D121.112   D121.221 D121.222   D121.121 D121.122
%   D122.111 D122.112   D122.221 D122.222   D122.121 D122.122
%
% INPUT
%  - U:       canonical displacement
%  - Ut:      fluctuation (also called localisator) (ut_ij)
%  - Utt:     localisator of 2nd order (utt_ijk)
%  - Ch0:     homogenized elasticity tensor Chom
%  - mesh:    pdetool mesh struct
%  - matprop: struct associated to the material properties
%  - psi:     level-set function
%  - volRVE:  total volume of the RVE
%
% OUTPUT
%  - Dh:      homogenized tensor Dhom
%
% HISTORY
% A.A. Novotny, V. Calisti  03/2020: code implementation.
%**************************************************************************

function Dh = tens_Dh(U,Ut,Utt,Ch0,mesh,matprop,psi,volRVE)

gamma = matprop.gamma; la0 = matprop.la0; mu0 = matprop.mu0;
% element characteristic function: 1 if psi > 0 and gama if psi <= 0
p  = mesh.p;    t = mesh.t;
% for the scalar product of matrix fields in Voigt writing
vcof = [1 1 2] ; C0 = matprop.C0 ;
es1 = matprop.es1 ; es2 = matprop.es2 ; 

tgamma = pdeintrp(p,t,(psi<0) + gamma*(psi>=0));
aux  = tgamma.*mesh.area;

% Calcul of Dhom 
z = 1; %compteur 

for i=1:3
    % i = ij  (1=11,2=22,3=12)
    
    % utt_ijk
    %   i1 = ij(k=1)
    Utti1  = pdeintrp(p,t,Utt(:,2*i-1));
    %   i2 = ij(k=2)
    Utti2  = pdeintrp(p,t,Utt(:,2*i  ));
    
    % sigma(utt_ij(k=1))
    stti1 = getstress(Utt(:,2*i-1),p,t,la0,mu0);
    % sigma(utt_ij(k=2))
    stti2 = getstress(Utt(:,2*i  ),p,t,la0,mu0);
    
    % ut_ij
    Uti  = pdeintrp(p,t,Ut(:,i));
    
    %~~~~
    % sigma(u_ij)
    si   = getstress(U(:,i),p,t,la0,mu0);
    %~~~~
    
    for j=1:3
        % j = pq  (1=11,2=22,3=12)
        
        % sigma(u_pq)
        sj   = getstress(U(:,j),p,t,la0,mu0);
        
        % ut_pq
        Utj  = pdeintrp(p,t,Ut(:,j));
        
        %~~~~
        % utt_pqr
        %   j1 = pq(r=1)
        Uttj1  = pdeintrp(p,t,Utt(:,2*j-1));
        %   j2 = pq(r=2)
        Uttj2  = pdeintrp(p,t,Utt(:,2*j  ));
        
        % eps(utt_pq(r=1))
        ettj1 = getstrain(Utt(:,2*j-1),p,t );
        % eps(utt_pq(r=2))
        ettj2 = getstrain(Utt(:,2*j  ),p,t );
        %~~~~
        
        %**********************************************************
        % Expression SYMMETRIC
        % (k,r)=(1,1)----------------------------------------------
        d = aux.*( vcof*((stti1 + C0*es1*Uti).*(ettj1 + es1*Utj)) ...
                  - vcof*(si.*(es1*Uttj1))...
                  - vcof*(sj.*(es1*Utti1)) );
        dh(z) = sum(d,2);
        z=z+1;
        % (k,r)=(2,2)----------------------------------------------
        d = aux.*( vcof*((stti2 + C0*es2*Uti).*(ettj2 + es2*Utj)) ...
                  - vcof*(si.*(es2*Uttj2))...
                  - vcof*(sj.*(es2*Utti2)) );
        dh(z) = sum(d,2);
        z=z+1;
        % (k,r)=(1,2)----------------------------------------------
        d = aux.*( vcof*((stti1 + C0*es1*Uti).*(ettj2 + es2*Utj)) ...
                  - vcof*(si.*(es1*Uttj2))...
                  - vcof*(sj.*(es2*Utti1)) );
        dh(z) = sum(d,2);
        z=z+1;
        % (k,r)=(2,1)----------------------------------------------
        d = aux.*( vcof*((stti2 + C0*es2*Uti).*(ettj1 + es1*Utj)) ...
                  - vcof*(si.*(es2*Uttj1))...
                  - vcof*(sj.*(es1*Utti2)) );
        dh(z) = sum(d,2);
        z=z+1;

      
    end
  
    
end

dh  = (1/volRVE)*dh;

% Write Dhom in matrix form :
%   D111.111 D111.112   D111.221 D111.222   D111.121 D111.122
%   D112.111 D112.112   D112.221 D112.222   D112.121 D112.122
%
%   D221.111 D221.112   D221.221 D221.222   D221.121 D221.122
%   D222.111 D222.112   D222.221 D222.222   D222.121 D222.122
%
%   D121.111 D121.112   D121.221 D121.222   D121.121 D121.122
%   D122.111 D122.112   D122.221 D122.222   D122.121 D122.122
D = zeros(6,6);
t = 1;
for i = 1:3
    for j = 1:3
        D(2*i-1,2*j-1) = dh(t) ;
        t = t + 1;
        D(2*i  ,2*j  ) = dh(t) ;
        t = t + 1;
        D(2*i-1,2*j  ) = dh(t) ;
        t = t + 1;
        D(2*i  ,2*j-1) = dh(t) ;
        t = t + 1;
    end
end

Dh.tens = D;
% Dh.eig  = eig(D);  
Dh.nm   = norm(D);

end
