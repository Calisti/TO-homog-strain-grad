%***********************************************************************
% Topological Derivative computation associated to the homogenised
% tensor D^hom
%***********************************************************************
% INPUT
%  - mesh:    pdetool mesh struct
%  - U:       canonical displacement
%  - Ut:      fluctuation (also called localisator) (ut_ij)
%  - Utt:     localisator of 2nd order (utt_ijk)
%  - P:       adjoint solutions p_ijk^r for calculus of DtD
%  - matprop: struct associated to the material properties
%  - psi:     level-set function
%  - volRVE:  total volume of the RVE
%
% OUTPUT
%  - DtD,DtDm:  topological derivative of D^hom, the first one is
%               written as a vector field, the second one is written as
%               a matrix field
%
% HISTORY
% A.A. Novotny, V. Calisti  03/2020: code implementation.
%***********************************************************************

function [DtD,DtDm] = topoD(mesh,U,Ut,Utt,Ch,K,bc,matprop,psi,cost,volRVE)

% loading the adjoint sates
P = adjD(Ut, Utt, K, psi, mesh, matprop, bc,volRVE);
% cost properties
IDv = cost.Idx.IDv;
% mesh properties
p = mesh.p; t = mesh.t; nt = size(t,2);
tarea  = mesh.area; 
% material properties
E0 = matprop.E0; gamma = matprop.gamma; 
la0 = matprop.la0; mu0 = matprop.mu0;
% help for computing
tr  = [1 1 0]; %n
vc  = [1 1 2] ; C0 = matprop.C0 ; %n
es1 = matprop.es1 ; es2 = matprop.es2 ; %n 
ES  = [es1 ; es2]; 

% element characteristic function: 1 if psi > 0 and gama if psi <= 0
tgamma       = pdeintrp(p,t,(psi<0) + gamma*(psi>=0));
tgammasquare = tgamma.*tgamma ;
tarea3       = [tgamma.*tarea ; tgamma.*tarea ; tgamma.*tarea];

tchi = pdeintrp(p,t,(psi < 0));
measchi = sum((tarea).*tchi);  coef_chi = 0 ;

% characteristic function for smoothing topo der
pchi = (psi<0); tchi = pdeintrp(p,t,pchi);

% extra parameters for topological derivative calculation
tE = E0*tgamma; alpha = (la0+mu0)/mu0; beta = (la0+3*mu0)/(la0+mu0);

% 1. topological derivative at the bulk fase 
cf1 = - ((1-gamma)/(1+beta*gamma))./tE; 
cf2 = - cf1.*( ( 1 - gamma*(beta-2*alpha) )/(1+alpha*gamma ) ); 
cf1 = cf1.*tgammasquare ;
cf2 = cf2.*tgammasquare ;
cf3 = (1/measchi)*coef_chi ; 
% 2. topological derivative at the inclusion
gammai = 1/gamma;
cf1i = - ((1-gammai)/(1+beta*gammai))./tE; 
cf2i = - cf1i.*( ( 1 - gammai*(beta-2*alpha) )/(1+alpha*gammai ) ); 
cf1i = cf1i.*tgammasquare ;
cf2i = cf2i.*tgammasquare ;
cf3i = (1/measchi)*coef_chi ; 

% Init Loop
nIv = length(IDv);
Dte = zeros(nIv,nt);   % top.der. : at bulk fase
ze   = 1 ; % counter   
Dti = zeros(nIv,nt);   % top.der. : at inclusion
zi   = 1 ; % counter 
Dts = zeros(nIv,nt);   % top.der. : smoothed and signed
zs   = 1 ; % counter 
% Define (k,r)
KR = [ 1 , 1 , 1 ;...
       2 , 2 , 2 ;...
       3 , 1 , 2 ;...
       4 , 2 , 1     ];
RK = [ 1 , 2 , 4 , 3 ];

for count = 1:nIv
    
    Nl = IDv(count);
    % i = ij
    i  = floor( (Nl-1)/12 ) + 1;
    % j = pq
    j  = rem(  (floor( (Nl-1)/4  ) + 1) - 1, 3) + 1;
    % kr 
    Nl4 = rem(Nl-1,4)+1; k = KR( Nl4, 2); r = KR( Nl4, 3);
    esk = ES( 3*(k-1)+1 : 3*(k-1)+3 , : );
    esr = ES( 3*(r-1)+1 : 3*(r-1)+3 , : );
    
    % i = ij
    
    %***utt_ijk
    Uttik = pdeintrp(p,t,Utt(:,2*i + (k-2)));
    %***sigma(utt_ijk)
    sttik = getstress(Utt(:,2*i + (k-2)),p,t,la0,mu0);
    %***ut_ij
    Uti  = pdeintrp(p,t,Ut(:,i));
    %***sigma(u_ij)
    si   = getstress(U(:,i),p,t,la0,mu0);
    %***sigma(p_ijk^r)
    spikr = getstress(P(:,4*i + (Nl4-4)),p,t,la0,mu0);
    
    % j =pq
    
    %***utt_pqr
    Uttjr = pdeintrp(p,t,Utt(:,2*j + (r-2)));
    %***sigma(utt_pqr)
    sttjr = getstress(Utt(:,2*j + (r-2)),p,t,la0,mu0);
    %***ut_pq
    Utj  = pdeintrp(p,t,Ut(:,j));
    %***sigma(u_pq)
    sj   = getstress(U(:,j),p,t,la0,mu0);
    %***sigma(p_pqr^k)
    spjrk=getstress(P(:,4*j + (RK(Nl4)-4)),p,t,la0,mu0);

    
    %******************************************************************
    % CALCULI coeff D^hom
    
    % 1. BULK PHASE ~~~~~~~~~~~~~~~~
    dte = 4*cf1.*( vc*(- si.*(spjrk + C0*esk*Uttjr)                   ...
                       - sj.*(spikr + C0*esr*Uttik)                   ...
                       +(sttik + C0*esk*Uti).*(sttjr + C0*esr*Utj) ) )...
          + cf2.*(-(tr*si).*(tr*(spjrk + C0*esk*Uttjr))               ...
                  -(tr*sj).*(tr*(spikr + C0*esr*Uttik))               ...
                  +(tr*(sttik+C0*esk*Uti)).*(tr*(sttjr+C0*esr*Utj)) ) ...
          + cf3 *( vc*(-volRVE*diag(Ch(j,:))*(esr*Uttik)              ...
                       -volRVE*diag(Ch(i,:))*(esk*Uttjr)              ...
                +diag(sum((tarea3).*(sttjr+C0*esr*Utj),2))*(esk*Uti)  ...
                +diag(sum((tarea3).*(sttik+C0*esk*Uti),2))*(esr*Utj) ) );

    
    Dte(ze,:) = dte;  ze = ze + 1 ;
    
    % 2. INCLUSION ~~~~~~~~~~~~~~~~
    dti = 4*cf1i.*( vc*(- si.*(spjrk + C0*esk*Uttjr)                   ...
                        - sj.*(spikr + C0*esr*Uttik)                   ...
                        +(sttik + C0*esk*Uti).*(sttjr + C0*esr*Utj) ) )...
          + cf2i.*(-(tr*si).*(tr*(spjrk + C0*esk*Uttjr))               ...
                   -(tr*sj).*(tr*(spikr + C0*esr*Uttik))               ...
                   +(tr*(sttik+C0*esk*Uti)).*(tr*(sttjr+C0*esr*Utj)) ) ...
          + cf3i *( vc*(-volRVE*diag(Ch(j,:))*(esr*Uttik)              ...
                        -volRVE*diag(Ch(i,:))*(esk*Uttjr)              ...
                 +diag(sum((tarea3).*(sttjr+C0*esr*Utj),2))*(esk*Uti)  ...
                 +diag(sum((tarea3).*(sttik+C0*esk*Uti),2))*(esr*Utj) ) );
    
    Dti(zi,:) = dti;  zi = zi + 1 ;
    
    % 3. SMOOTHING + SIGN ~~~~~~~~~
    Dts(zs,:) = - tchi.*dte + (1-tchi).*dti;
    zs = zs + 1 ;       
    
end

Dts = Dts/volRVE ;
DtDm = [];

%--Only saving the needed values
DtD = Dts ;

end

