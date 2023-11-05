%**************************************************************************
% Topological Derivative a Functional with Volume Constraint  
%**************************************************************************
% DESCRIPTION
% Computes the topological derivative
% 
% INPUT
%  - mesh:     pdetool mesh struct
%  - U:        pdetool solution
%  - Ut:       first order corrector
%  - Utt:      second order correcor
%  - K:        stiffness matrix of the reduced system (periodic bc)
%  - bc:       extra boundry conditions struct
%  - volume:   current volume
%  - matprop:  material properties struct
%  - psi:      level-set function
%  - params:   topology optimization parameters struct
%  - pdecoef:  pdetool coeficiets struct
%  - volRVE:   total volume of the RVE
%  - cost:     option for the choice of the shape functional
%  - Ch,Dh:    homogenized tensors
%
% OUTPUT
%  - dt: topological derivative -> dt = -DT (bulk) and dt = DT (inclusion)
%
% HISTORY
% S. Amstutz & A.A. Novotny      06/2009: code implementation.
% D.E. Campe�o                   12/2010: code updating.
% J-M.C. Farias
% A.A. Novotny                   06/2012: code updating.
% A.A. Novotny, V. Calisti          2020: code updating.
%**************************************************************************

function dt = topder(mesh,U,Ut,Utt,K,bc,volume,matprop,psi,params, ...
                                                    volRVE,cost,Ch,Dh)

% INITIALISATION
% #########################################################################

% mesh functions
p  = mesh.p; t = mesh.t;
nt = size(t,2);
% loading Ch and Dh
Ch = Ch.tens ;
Dh = Dh.tens ;
% loading top der of Ch
Dt = topoC(mesh,U,matprop,psi,cost,volRVE);

% load penalization parameters and init value of the functional for
% nomalization
penalization = params.penalization; penalty = params.penalty;
volfrac = params.volfrac; auglag = params.auglag; volinit = params.volinit;
sfJ0 = params.sfJ0; voltarget = volinit*volfrac;

% Compute the topological derivatives of coeff according to chosen costJ
% Load the top der of Dh
%  + eventually compute the topological derivative of Sh = Ch^{-1}
[DtDhom,~] = topoD(mesh,U,Ut,Utt,Ch,K,bc,matprop,psi,cost,volRVE);
DtCm       = Dt.DtCm;
IS = cost.Idx.IS ;
Sh = 0;
if ~ isempty(IS)
    Sh     = inv(Ch);
    MinvCh = [ Sh(1,1)*speye(nt), Sh(1,2)*speye(nt), Sh(1,3)*speye(nt);...
        Sh(2,1)*speye(nt), Sh(2,2)*speye(nt), Sh(2,3)*speye(nt);...
        Sh(3,1)*speye(nt), Sh(3,2)*speye(nt), Sh(3,3)*speye(nt)   ];
    DtSm = - Sh*( DtCm*MinvCh ) ;
end   

% COMPUTE THE TOPOLOGICAL DERIVATIVE
% #########################################################################
      
GradFunJ = cost.GradFunJ;
Idx      = cost.Idx;
lC   = Idx.lC;   lS   = Idx.lS;
IC   = Idx.IC;   IS   = Idx.IS;   IDm  = Idx.IDm;

D1 = [] ; D2 = [] ;
nC = size(lC,1) ; nS = size(lS,1);
for i=1:nC
    D1(i,:) = DtCm( lC(i,1), 1+ (lC(i,2)-1)*nt : (lC(i,2))*nt ) ;
end
for i=1:nS
    D2(i,:) = DtSm( lS(i,1), 1+ (lS(i,2)-1)*nt : (lS(i,2))*nt ) ;
end

DTens = [ D1 ; D2 ; DtDhom ];

dt = sum( ( GradFunJ( [ Ch(IC); Sh(IS); Dh(IDm) ] ) ).*DTens , 1) ;
dt = pdeprtni(p,t,dt);
  

% #########################################################################
% Include the volume constraint sensitivity

dt = dt / sfJ0;
if penalization == 1
    dt = dt + penalty / volinit;
    %         figure(20); pdesurf(p,t,dt);
    %         pause
elseif penalization == 2
    coef = volume / voltarget;
    dt = dt + penalty*((coef <= 1)*0.0 - (coef > 1)*2.0*(1.0-coef)) / voltarget;
elseif penalization == 3
    coef = volume / voltarget;
    dt = dt + max(0,penalty - auglag*(1.0-coef)) / voltarget;
end

end
