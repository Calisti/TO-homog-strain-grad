%**************************************************************************
% Solution update   
%**************************************************************************
% DESCRIPTION
% Updates solution of the multi-scale pde linear system. 
% 
% INPUT
%  - Ur:   pdetool reduced solution
%  - Us:   excitation of multi-scale system
%  - bc:   extra boundary conditions struct
%
% OUTPUT
%  - U:    updated pdetool solution (and pre-load vector)
%
% HISTORY
% S.M. Giusti         02/2011: code implementation.
% A.A. Novotny        06/2012: code updating.
%**************************************************************************

function U = solupdate(Ur,Us,bc)

if not(size(Ur,2)==size(Us,2))
    error('First two inputs of solupdate do not have the same size')
end
aux = size(Ur,2);

nbdofa = length(bc.dofs.all);
nbdofi = length(bc.dofs.i);
nbdofp = length(bc.dofs.p);
nbdofm = length(bc.dofs.m);
U = zeros(nbdofa,aux);
U(bc.dofs.i,:) = Ur(1:nbdofi,:) + Us(1:nbdofi,:);
U(bc.dofs.p,:) = Ur(nbdofi+1:length(Ur),:) + Us(nbdofi+1:nbdofi+nbdofp,:);
U(bc.dofs.m,:) = Ur(nbdofi+1:length(Ur),:) + Us(nbdofi+nbdofp+1:nbdofi+nbdofp+nbdofm,:);
U(bc.dofs.c,:) = Us(nbdofi+nbdofp+nbdofm+1:nbdofa,:);
    
end
