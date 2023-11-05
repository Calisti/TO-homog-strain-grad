%**************************************************************************
% Derivative matrix 
%**************************************************************************
% DESCRIPTION
% np = nbr of nodes , nt = nbr of triangles
% for z = x,y  
% dz = (nt,np)matrix such that f=(np,1)vector field defined on the nodes : 
%   dz'*f = (df)/(dz) = (1,nt)vector 
%
% CALLED in ...
% ... assemFload.m
%
% INPUT
%  - mesh:     pdetool mesh struct
%
% OUTPUT
%  - [dx,dy] = matrices defined in description
%
% HISTORY
% A.A. Novotny, V. Calisti : 03/2020 
%**************************************************************************

function [dx,dy] = sfderivatives(mesh)
   
    p = mesh.p; t = mesh.t; np = size(p,2); nt = size(t,2);    
    [~,e1x,e1y,e2x,e2y,e3x,e3y]=pdetrg(p,t);
    
    k = 1:nt;
    dx = sparse(k,t(1,k),e1x,nt,np);
    dx = dx+sparse(k,t(2,k),e2x,nt,np);
    dx = dx+sparse(k,t(3,k),e3x,nt,np);
    dy = sparse(k,t(1,k),e1y,nt,np);
    dy = dy+sparse(k,t(2,k),e2y,nt,np);
    dy = dy+sparse(k,t(3,k),e3y,nt,np);
    
    %dsf = [dx,dy];
    
end