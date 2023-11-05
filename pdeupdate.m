%**************************************************************************
% Pdetool update   
%**************************************************************************
% DESCRIPTION
% Updates pdetool linear system, written for the complete geometry of the
% unit cell, into the degrees of freedom used for periodic problem.
% 
% CALLED in :
% ...solve.m
%
% INPUT
%  - K:    pdetool stiffness matrix
%  - bc:   extra boundary conditions struct
%  - mesh: mesh parameters struct 
%
% OUTPUT
%  - [Kr,Fr,US]: updated pdetool stiffness matrix, load vector,
%                and linear part of the micro displacement
%
% HISTORY
% J-M.C. Farias             12/2010: code implementation.
% D.E. Campeao     
% A.A. Novotny
%**************************************************************************

function [Kr,Fr,Us] = pdeupdate(K,bc,mesh)

dofs = bc.dofs;
nodes = mesh.nodes;
p = mesh.p; t = mesh.t;

% Set relevant stiffness sub-matrices
Kii = K(dofs.i(:),dofs.i(:));
Kip = K(dofs.i(:),dofs.p(:));
Kim = K(dofs.i(:),dofs.m(:));
Kic = K(dofs.i(:),dofs.c(:));
Kpi = K(dofs.p(:),dofs.i(:));
Kpp = K(dofs.p(:),dofs.p(:));
Kpm = K(dofs.p(:),dofs.m(:));
Kpc = K(dofs.p(:),dofs.c(:));
Kmi = K(dofs.m(:),dofs.i(:));
Kmp = K(dofs.m(:),dofs.p(:));
Kmm = K(dofs.m(:),dofs.m(:));
Kmc = K(dofs.m(:),dofs.c(:));

% Assemble reduced system forcing term
Kstar = [    Kii      Kip      Kim      Kic     ;
    Kpi+Kmi  Kpp+Kmp  Kpm+Kmm  Kpc+Kmc  ];

global_basis = eye(3);
Fr=[]; Us=[];
for i=1:3
    strain_tensor = [  global_basis(1,i)     global_basis(3,i)/2 ;
        global_basis(3,i)/2   global_basis(2,i)  ];

    Usxi = strain_tensor(1,1)*p(1,nodes.i)' + strain_tensor(1,2)*p(2,nodes.i)';
    Usyi = strain_tensor(2,1)*p(1,nodes.i)' + strain_tensor(2,2)*p(2,nodes.i)';
    Usxp = strain_tensor(1,1)*p(1,nodes.p)' + strain_tensor(1,2)*p(2,nodes.p)';
    Usyp = strain_tensor(2,1)*p(1,nodes.p)' + strain_tensor(2,2)*p(2,nodes.p)';
    Usxm = strain_tensor(1,1)*p(1,nodes.m)' + strain_tensor(1,2)*p(2,nodes.m)';
    Usym = strain_tensor(2,1)*p(1,nodes.m)' + strain_tensor(2,2)*p(2,nodes.m)';
    Usxc = strain_tensor(1,1)*p(1,nodes.c)' + strain_tensor(1,2)*p(2,nodes.c)';
    Usyc = strain_tensor(2,1)*p(1,nodes.c)' + strain_tensor(2,2)*p(2,nodes.c)';

    Usi = [Usxi;Usyi;Usxp;Usyp;Usxm;Usym;Usxc;Usyc];
    clear Usxi Usyi Usxp Usyp Usxm Usym Usxc Usyc;
    Fri = - Kstar * Usi;
    Fr  = [Fr,Fri]; Us = [Us,Usi];

end

clear Kstar Fri Usi;
% Assemble reduced system matrix
Kr = [   Kii         Kip+Kim      ;
    Kpi+Kmi   Kpp+Kpm+Kmp+Kmm ];


end
