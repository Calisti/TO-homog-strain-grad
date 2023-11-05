%**************************************************************************
% Compute the first order correctors  
%**************************************************************************
% DESCRIPTION
% Assemble and solve the pdetool finite element linear system 
% for the equations of ut_ij = FIRST ORDER CORRECTOR,
% (saved in Ut(size=(3,np))), and with 
% u_ij = ut_ij + C*(ei x_s ej) (saved in U(size=(3,np)))
% Ut = [ Ut11 Ut22 Ut12 ]  (=Ut_ij) 
% U  = [  U11  U22  U12 ]  (=U_ij) 
% 
% INPUT
%  - psi:     level-set function
%  - mesh:    pdetool mesh struct
%  - matprop: material properties struct
%  - pdecoef: pdetool coeficients struct
%  - bc:      extra boundry conditions struct
%
% OUTPUT
%  - [U,Ut,K]: micro displacement, fluctuation (also called localisator)
%              and stiffness matrix of the reduced system (periodic bc)  
%
% HISTORY
% J-M.C. Farias             12/2010: code implementation.
% A.A. Novotny              12/2010: code updating.
% S.M. Giusti               02/2011: code updating.
% A.A. Novotny, V. Calisti  03/2020: code updating.
%**************************************************************************

function [U,Ut,K] = solve(psi, mesh, matprop, pdecoef, bc)
 
    gamma  = matprop.gamma; 
    % mesh functions and pde coeff
    p = mesh.p; t = mesh.t; a = pdecoef.a; f = pdecoef.f; c0 =  pdecoef.c;
    np = size(p,2);
    
    % element characteristic function: 1 if psi > 0 and gama if psi <= 0
    tgamma = pdeintrp(p,t,(psi<0) + gamma*(psi>=0));
    
    c = c0 * tgamma; % effective elasticity tensor 
    
    [K,~,~] = assema(p,t,c,a,f);
    % modify system to solve the pbm depending on bc
    [K,Fload,Us] = pdeupdate(K,bc,mesh);
    % solve linear system
    Ut = K \ Fload;
    % extend solution [ Ui  U+ ] to all nodes : [ Ui  U+  U+  0 ]
    % with nodes : [ interior  faces+  faces-  corners ]
    U = solupdate(Ut,Us,bc);
    Ut = solupdate(Ut,zeros(2*np,3),bc);
    % fix the average of Ut equal to zero; 
    Ut = fixmean(Ut,psi,mesh,matprop,bc);  
    
end

