%**************************************************************************
% Local remesh
%**************************************************************************
% DESCRIPTION
% Refine the mesh only in the areas where the level-set value are high, 
% also refine it on the edges to not break the periodicity of the mesh
%
% INPUT
%  - mesh:       current mesh                
%  - psi & psi0: current & init. level-sets 
%
% OUTPUT
%  - p,e,t:      updated mesh                               
%  - psi & psi0: current & init. level-sets on the new mesh
%
% HISTORY
% V. Calisti  04/2021: code implementation.
%***********************************************************************

function [p,e,t,psi,psi0] = remeshlocal(mesh,psi,psi0)

% Load mesh 
g = mesh.g; p = mesh.p; e = mesh.e; t = mesh.t;

% On raffine dans une région incluant la zone des valeurs 
% de psi en dessous d'un certain seuil + raffine les bords
for i = 1:1
    tpsi = pdeintrp(p,t,psi);
    % Tentative : remailler l'interface
    % --Find triangles with sign change of psi -> e0
    Tpsi = sort( psi(t(1:3,:)) , 1);
    Tpsi = Tpsi(1,:).*Tpsi(3,:);
    e0   = find( Tpsi < 0 );
    % --To remesh dans un voisinage des valeurs minimales de psi
    minpsi = min(tpsi);
    % TWo gpes of triangles with rpct to a threshold 
    % el0 -> remesh,  elc -> no remesh
    el0 = find( tpsi <  ( (1/4)  )*minpsi );
    elc = find( tpsi >= ( (1/4)  )*minpsi );
    % Get the corners nodes of these triangles
    it10=t(1,el0); it20=t(2,el0); it30=t(3,el0);
    it1c=t(1,elc); it2c=t(2,elc); it3c=t(3,elc);
    % Get the xy coordinates of centers of triangles
    xt0 = (p(1,it10)+p(1,it20)+p(1,it30))/3; 
    yt0 = (p(2,it10)+p(2,it20)+p(2,it30))/3;  
    xtc = (p(1,it1c)+p(1,it2c)+p(1,it3c))/3; xtc = xtc';
    ytc = (p(2,it1c)+p(2,it2c)+p(2,it3c))/3; ytc = ytc'; 
    % Get the elc triangles in the vicinity of the el0 tr. for remeshing
    Mat = (xt0-xtc).^2 + (yt0-ytc).^2 < (0.05)^2 ; 
    ielc = find( sum(double(Mat),2) > 0 );
    el = [el0 , elc(ielc)]'; 
    % --To remesh the boundary
    PB = mesh.nodes.b';
    tb = sum(PB == t(1,:), 1) + sum(PB == t(2,:), 1) + sum(PB == t(3,:), 1);
    eb = find(tb == 1); 
    el = unique([el; eb'; e0']);
    % --Remesh
    [p,e,t,psi] = refinemesh(g,p,e,t,[psi psi0],el); 
    psi0 = psi(:,2); psi = psi(:,1);
end

end

