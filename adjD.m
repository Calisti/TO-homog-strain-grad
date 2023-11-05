%**************************************************************************
% Compute Adjoint State for Dh tensor topo.der. 
%**************************************************************************
% DESCRIPTION
% Assemble and solve the pdetool finite element linear system 
% for the equations of p_ijk^r = ADJOINT STATE for Dh
% These field is saved in P (size=(12,np)) , where np = nbr of nodes          
% P   = (Pijkr)  = [ P1111 P1122 P1112 P1121 , 
%                    P2211 P2222 P2212 P2221 ,
%                    P1211 P1222 P1212 P1221  ]
% 
% INPUT
%  - Ut:      first order corrector
%  - Utt:     second order corrector
%  - K:       stiffness matrix of the reduced system (periodic bc)
%  - psi:     level-set function
%  - mesh:    pdetool mesh struct
%  - matprop: material properties struct
%  - bc:      extra boundry conditions struct
%  - volRVE:  volume of the RVE
%
% OUTPUT
%  - P: pdetool solutions p_ijk^r
%
% HISTORY
% A.A. Novotny, V. Calisti  03/2020: code implementation.
%**************************************************************************

function P = adjD(Ut, Utt, K, psi, mesh, matprop, bc,volRVE)

gamma = matprop.gamma; la0 = matprop.la0; mu0 = matprop.mu0;
% usefull coefficient, help for computing
C0 = matprop.C0;
e1 = matprop.e1;     e2 = matprop.e2;
es1 = matprop.es1;   es2 = matprop.es2;
% mesh functions and pde coeff
p = mesh.p; t = mesh.t;
%a = pdecoef.a; c0 =  pdecoef.c;%old
tarea = [mesh.area;mesh.area];
% element characteristic function: 1 if psi > 0 and gama if psi <= 0
tgamma  = pdeintrp(p,t,(psi<0) + gamma*(psi>=0));
tgamma2 = [tgamma;tgamma]; tgamma3 = [tgamma;tgamma;tgamma];
% for the constant term of equation
np = size(p,2);
ni = size(mesh.nodes.i,2); nplus = size(mesh.nodes.p,2);
% normalized characteristic function
tchi = ones(1,size(t,2));
tchi2_nzd = [tchi ; tchi];

% -------------------------
% CALCULUS OF  p_ijk^r     represented by:    P
% -------------------------

F = zeros(2*(ni+nplus),12);

for i=1:3
    % i = ij   (1=11,2=22,3=12)
    
    % ut_ij
    Uti  = pdeintrp(p,t,Ut(:,i));
    % Fkr, for k=1,2 r=1,2
    [F11,F12] = assemFload(mesh,p,t,tgamma2,la0,mu0,Utt(:,2*i-1));
    [F21,F22] = assemFload(mesh,p,t,tgamma2,la0,mu0,Utt(:,2*i  ));
    % stress of utt_ijk :   k=1 s1   ,   k=2 s2
    si1 = getstress(Utt(:,2*i-1),p,t,la0,mu0);
    si2 = getstress(Utt(:,2*i  ),p,t,la0,mu0);
    
    % (k,r)=(1,1)--------------------------------------------
    f = tgamma3(1:2,:).*( e1*( C0*es1*Uti + si1 ) );
    moy = (1/volRVE)*sum(tarea.*f,2) ;
    f = f - diag(moy)*tchi2_nzd ;
    [~,~,Fi] = assema(p,t,[],[],f);
    Fi = Fi - F11;
    % reduce Fload for solving Rr
    Fi = reducedF(Fi,bc);
    F(:,4*i-3) = Fi;
    
    % (k,r)=(2,2)--------------------------------------------
    f = tgamma3(1:2,:).*( e2*( C0*es2*Uti + si2 ) );
    moy = (1/volRVE)*sum(tarea.*f,2) ;
    f = f - diag(moy)*tchi2_nzd ;
    [~,~,Fi] = assema(p,t,[],[],f);
    Fi = Fi - F22;
    % reduce Fload for solving Rr
    Fi = reducedF(Fi,bc);
    F(:,4*i-2) = Fi;
    
    % (k,r)=(1,2)--------------------------------------------
    f = tgamma3(1:2,:).*( e2*( C0*es1*Uti + si1 ) );
    moy = (1/volRVE)*sum(tarea.*f,2) ;
    f = f - diag(moy)*tchi2_nzd ;
    [~,~,Fi] = assema(p,t,[],[],f);
    Fi = Fi - F12;
    % reduce Fload for solving Rr
    Fi = reducedF(Fi,bc);
    F(:,4*i-1) = Fi;
    
    % (k,r)=(2,1)--------------------------------------------
    f = tgamma3(1:2,:).*( e1*( C0*es2*Uti + si2 ) );
    f = f - diag(moy)*tchi2_nzd ;
    [~,~,Fi] = assema(p,t,[],[],f);
    Fi = Fi - F21;
    % reduce Fload for solving Rr
    Fi = reducedF(Fi,bc);
    F(:,4*i  ) = Fi;
    
end

% solve the reduced problem
Pr = K \ F;
% extend solution Pr to P
P = zeros(2*np,12);
P = solupdate(Pr,zeros(2*np,size(Pr,2)),bc);
% impose P with mean value = 0
P = fixmean(P,psi,mesh,matprop,bc);
    
end

