%**************************************************************************
% Initialization of mscale.m code 
%**************************************************************************
% HISTORY
% S. Amstutz & A.A. Novotny      06/2009: code implementation.
% S.M. Giusti                    03/2011: code updating.
% A.A. Novotny                   06/2012: code updating.
% V.   Calisti                   2020-21: code updating.
%**************************************************************************

function [mesh, pdecoef, matprop, params, bc, psi0, cost] = init_file 

%%INPUT PARAMETERS
% rectangle unit cell dimensions
a = 1; b = 1; 
% initial mesh size
ni = 40; 
% s =  string expression of the functional, is written with the coeff
%      of the selected tensors *ij,.... . 
%      These coeff. have to be surrounded with empty spaces.
%      Also, tensors need to be written Thij with i <= j 
%      it is not restrictive because all tensors are symetric
%      ex : 
%            s = ' Ch11 * Sh23 / ( Dh44 ^2 ) ';
s = '- Dh22 / Ch11 '
% material properties: Young modulus, material contrast, poisson ratio
E0 = 1.0; gamma = 0.01 ; nu = 0.3;

%% geometry

p = [ 0 , 0;
      a , 0;
      a , b;
      0 , b];
g = [2; size(p,1); p(:,1); p(:,2)];

% composing the geometry
g = decsg(g);

%% algorithm parameters

% minimum allowed 'k'. Used by the line-search procedure
params.kmin = 1.0E-3;
% starting value of 'k' ( k <= 1 )
params.kstart = 1 ;
% stop criterion
params.stop = 1.0*pi/180;

% method of penalization
% (1) linear
% (2) exact quadratic
% (3) augmented-lagrangian
params.penalization = 1;

% penalty parameter   
params.penalty = 0 ;    % method 1, 2 and 3                                %= 0.6225;
params.volfrac = 0.75;  % method 2 and 3
params.voleps  = 0.01;  % method 2 and 3
params.auglag  = 0.0;   % method 3
params.epsilon = 0.01;  % method 3

%% mesh generation

[p,e,t] = poimesh(g,a*ni,b*ni); [p,e,t] = refinemesh(g,p,e,t,'longest');
ghold  = 0; nsteps = 0; remesh  = 'regular';

for i=1:2*nsteps
    [p,e,t] = refinemesh(g,p,e,t,remesh);
end

area = pdetrg(p,t);

mesh.remesh = remesh;
mesh.p = p;
mesh.e = e;
mesh.t = t;
mesh.g = g;
mesh.area  = area;
mesh.ghold = ghold;
mesh.ni = ni;
mesh.a  = a;  mesh.b  = b;

%% multi-scale model

bc.model = 'periodic';

%% cost function

[FunJ, GradFunJ, Idx] = init_func(s);
cost.FunJ             = FunJ;
cost.GradFunJ         = GradFunJ;
cost.Idx              = Idx;
cost.name             = s;
cost.beta = 0.5; 
% normalization option : (0 = no) or (1 = yes):
cost.normzn = 0;

%% material properties

matprop.E0 = E0; matprop.gamma = gamma; matprop.nu = nu;

%% pde coeficients

la0 = nu*E0/((1+nu)*(1-2*nu)); mu0 = E0/(2*(1+nu)); % plane strain
la0 = 2*mu0*la0/(la0+2*mu0); % plane stress
matprop.la0 = la0; matprop.mu0 = mu0;

% For calculus with voigt notation
matprop.vcoef = [1 0 0; 0 1 0 ; 0 0 2];            
matprop.vc    = [1 1 2];                           
matprop.C0    = [2*mu0+la0   la0       0    ; ...  
                    la0    2*mu0+la0   0    ; ... 
                     0        0      2*mu0 ];  
matprop.C0ij  = [2*mu0+la0   la0       0    ; ... 
                    la0    2*mu0+la0   0    ; ... 
                     0        0       mu0  ];  
matprop.es1   = [1 0 ; 0 0 ; 0   0.5];      
matprop.es2   = [0 0 ; 0 1 ; 0.5   0];       
matprop.eS    = [matprop.es1 ; matprop.es2];  
matprop.e1    = [1 0 0 ; 0 0 1];   
matprop.e2    = [0 0 1 ; 0 1 0];     
matprop.tr    = [1 0 0];        

c = zeros(10,1);
c(1,:) = 2*mu0 + la0; c(3,:) = mu0;  c(5,:) = mu0;
c(6,:) = la0; c(8,:) = mu0; c(10,:) = 2*mu0 + la0;

pdecoef.c = c;
pdecoef.a = zeros(4,1);
pdecoef.f = zeros(2,1);

%% initial guess

% Definition of PSI0, the initial level-set function

%-------- Trou centré pour carré 1x1
aux = cos(pi*(p(1,:)-0.5)).^2 .* cos(pi*(p(2,:)-0.5)).^2;
psi0 = (aux-0.5)';


end
