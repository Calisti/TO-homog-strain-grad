%***********************************************************************
% Topological Derivative computation associated to the homogenised
% tensor Ch
%***********************************************************************
% DESCRIPTION
% Given the canonical stress tensor, computes the associated topological 
% derivative the homogenised elasticity tensor
% 
% INPUT
%  - p: pdetool points (nodes)
%  - t: pdetool triangles (elements)
%  - U: canonical displacement
%  - matprop: struct associated to the material properties
%  - psi: level-set function
%  - volRVE: total volume of the RVE
%
% OUTPUT
%  - Dt: topological derivative tensor field of Chom
%
% HISTORY
% S. Amstutz & A.A. Novotny      06/2009: code implementation.
%***********************************************************************

function Dt = topoC(mesh,U,matprop,psi,cost,volRVE)  

% load mesh
p = mesh.p; t = mesh.t;
% material properties
E0 = matprop.E0; gamma = matprop.gamma; 
la0 = matprop.la0; mu0 = matprop.mu0;
% help for comuting
vcoef = [1 1 2];
tr = [1 1 0];

% element characteristic function: 1 if psi > 0 and gama if psi <= 0
tgamma = pdeintrp(p,t,(psi<0) + gamma*(psi>=0));

% extra parameters for topological derivative calculation
tE = E0*tgamma; alpha = (la0+mu0)/mu0; beta = (la0+3*mu0)/(la0+mu0);

% nominal stress for each canonical direction 
s1 = getstress(U(:,1),p,t,la0,mu0); %[[1 0]; [0 0]]     : s(u_11)
s2 = getstress(U(:,2),p,t,la0,mu0); %[[0 0]; [0 1]]     : s(u_22)
s3 = getstress(U(:,3),p,t,la0,mu0); %[[0 1/2]; [1/2 0]] : s(u_12)

% effective stress
tgamma3=[tgamma;tgamma;tgamma]; 
s1=s1.*tgamma3; s2=s2.*tgamma3; s3=s3.*tgamma3;

% topological derivative at the bulk fase
coef1 = - ((1-gamma)/(1+beta*gamma))./tE; %new
coef2 = - coef1.*( ( 1 - gamma*(beta-2*alpha) )/(1+alpha*gamma) ); %new

dte1=coef1.*(4*vcoef*(s1.*s1)) + coef2.*(tr*s1).*(tr*s1); %Dt1111
dte2=coef1.*(4*vcoef*(s1.*s2)) + coef2.*(tr*s1).*(tr*s2); %Dt1122
dte3=coef1.*(4*vcoef*(s1.*s3)) + coef2.*(tr*s1).*(tr*s3); %Dt1112
dte4=coef1.*(4*vcoef*(s2.*s2)) + coef2.*(tr*s2).*(tr*s2); %Dt2222
dte5=coef1.*(4*vcoef*(s2.*s3)) + coef2.*(tr*s2).*(tr*s3); %Dt2212
dte6=coef1.*(4*vcoef*(s3.*s3)) + coef2.*(tr*s3).*(tr*s3); %Dt1212

% topological derivative at the inclusion
gamma = 1/gamma;
coef1 = - ((1-gamma)/(1+beta*gamma))./tE; %new
coef2 = - coef1.*( ( 1 - gamma*(beta-2*alpha) )/(1+alpha*gamma) ); %new

dti1=coef1.*(4*vcoef*(s1.*s1)) + coef2.*(tr*s1).*(tr*s1); %Dt1111
dti2=coef1.*(4*vcoef*(s1.*s2)) + coef2.*(tr*s1).*(tr*s2); %Dt1122
dti3=coef1.*(4*vcoef*(s1.*s3)) + coef2.*(tr*s1).*(tr*s3); %Dt1112
dti4=coef1.*(4*vcoef*(s2.*s2)) + coef2.*(tr*s2).*(tr*s2); %Dt2222
dti5=coef1.*(4*vcoef*(s2.*s3)) + coef2.*(tr*s2).*(tr*s3); %Dt2212
dti6=coef1.*(4*vcoef*(s3.*s3)) + coef2.*(tr*s3).*(tr*s3); %Dt1212

% smoothing of the topological derivative
% and sign to difine a target level-set
pchi = (psi<0); tchi = pdeintrp(p,t,pchi);

nt = size(t,2);
D = zeros(3,3*nt);
D(1,1     :nt  ) = (- tchi.*dte1 + (1-tchi).*dti1)/volRVE;% Dt1111;
D(2,1     :nt  ) = (- tchi.*dte2 + (1-tchi).*dti2)/volRVE;% Dt1122;
D(3,1     :nt  ) = (- tchi.*dte3 + (1-tchi).*dti3)/volRVE;% Dt1112;
D(1,nt+1  :2*nt) = (- tchi.*dte2 + (1-tchi).*dti2)/volRVE;% Dt1122;
D(2,nt+1  :2*nt) = (- tchi.*dte4 + (1-tchi).*dti4)/volRVE;% Dt2222;
D(3,nt+1  :2*nt) = (- tchi.*dte5 + (1-tchi).*dti5)/volRVE;% Dt2212;
D(1,2*nt+1:3*nt) = (- tchi.*dte3 + (1-tchi).*dti3)/volRVE;% Dt1112;
D(2,2*nt+1:3*nt) = (- tchi.*dte5 + (1-tchi).*dti5)/volRVE;% Dt2212;
D(3,2*nt+1:3*nt) = (- tchi.*dte6 + (1-tchi).*dti6)/volRVE;% Dt1212;
Dt.DtCm = D ;

end





