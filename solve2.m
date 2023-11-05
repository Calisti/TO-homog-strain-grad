%**************************************************************************
% Compute the second order correctors 
%**************************************************************************
% DESCRIPTION
% Assemble and solve the pdetool finite element linear system 
% for the equations of utt_ijk = SECOND ORDER CORRECTOR, 
% This field is respectively 
% saved in Utt (size=(6,np)) where np = nbr of nodes          
% Utt = (Uttijk) = [ Utt111 Utt112 Utt221 Utt222 Utt121 Utt122 ]
% 
% INPUT
%  - U, Ut:   micro displacement, and fluctuation (also called corrector)
%  - Ch:      homogenized elasticity tensor
%  - K:       stiffness matrix of the reduced system (periodic bc)
%  - psi:     level-set function
%  - mesh:    pdetool mesh struct
%  - matprop: material properties struct
%  - pdecoef: pdetool coeficients struct
%  - bc:      extra boundry conditions struct
%
% OUTPUT
%  - Utt:     pdetool solutions u_ijk
%
% HISTORY
% A.A. Novotny, V. Calisti  03/2020: code implementation.
%**************************************************************************

function Utt = solve2(U, Ut, Ch, K, psi, mesh, matprop, bc,volRVE)

gamma = matprop.gamma; la0 = matprop.la0; mu0 = matprop.mu0;
% usefull coefficient, help for computing
e1 = matprop.e1;     e2 = matprop.e2;
% mesh functions and pde coeff
p = mesh.p; t = mesh.t;
tarea = [mesh.area;mesh.area;mesh.area];%NEM
% element characteristic function: 1 if psi > 0 and gama if psi <= 0
tgamma  = pdeintrp(p,t,(psi<0) + gamma*(psi>=0));
tgamma2 = [tgamma;tgamma]; tgamma3 = [tgamma;tgamma;tgamma];
% effective elasticity tensor
Ch = Ch.tens ;
%c = c0 * tgamma;
% for the constant term of equation
%nt = size(t,2);%old
np = size(p,2);
ni = size(mesh.nodes.i,2); nplus = size(mesh.nodes.p,2);
% normalized characteristic function
tchi = ones(1,size(t,2));
tchi3_nzd = [tchi ; tchi; tchi];

% -------------------------------------
% CALCULUS OF  u_ijk                   represented by:   Utt
% -------------------------------------

% Calculus of the load Fload
Fload = zeros(2*(ni+nplus),6);

for i=1:3
    % i = ij  (1=11,2=22,3=12)
    
    s = getstress(U(:,i),p,t,la0,mu0); s = s.*tgamma3; 
    % Fk, for k=1,2
    [F1,F2] = assemFload(mesh,p,t,tgamma2,la0,mu0,Ut(:,i));
    % k = 1-----------------------------------------ek = (1,0)
    f = e1*( s - diag(Ch(:,i))*tchi3_nzd );
    [~,~,Fi] = assema(p,t,[],[],f);
    Fi = Fi - F1;
    % reduce Fload for solving Uttr
    Fi = reducedF(Fi,bc);
    Fload(:,2*i-1) = Fi;
    % k = 2-----------------------------------------ek = (0,1)
    f = e2*( s - diag(Ch(:,i))*tchi3_nzd );
    [~,~,Fi] = assema(p,t,[],[],f);
    Fi = Fi - F2;
    % reduce Fload for solving Uttr
    Fi = reducedF(Fi,bc);
    Fload(:,2*i) = Fi;
    
end

% solve the reduced problem
Uttr = K \ Fload;
% extend solution Uttr to Utt
Utt = zeros(2*np,6);
Utt = solupdate(Uttr,zeros(2*np,size(Uttr,2)),bc);
% impose mean value = 0
Utt = fixmean(Utt,psi,mesh,matprop,bc);

    
end

