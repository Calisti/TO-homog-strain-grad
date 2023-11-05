%**************************************************************************
% INITIALIZE the shape FUNCTIONAL
%**************************************************************************
% DESCRIPTION
% Take a string expression and calculate the associated cost functional 
% and its gradient, and calculate useful indices for the numerical   
% evaluation of the cost and its topological derivative.
%
% INPUT
%  - s = string expression of the functional, is written with the coeff
%        of the selected tensors *ij,.... . 
%        These coeff. have to be surrounded with empty spaces.
%        ex : 
%           s = ' Ch11 * Sh23 / ( Dh44 ^2 ) ';
%
% OUTPUT
%  - FunJ     = cost function (@ function handle)
%  - GradFunJ = gradient of the cost function (@ function handle)
%  - Idx      = indices (for identification of the selected coefficients
%               of homogenized tensors: * 
%               and their topological derivatives: Dt* )
%               Idx contains :
%               - arrays l*,... containing the indices of selected coeff
%                 of tensors *,...  written this way : * = (*)_ij
%                 example : if *11 and *12 are selected for the cost, then
%                           l* = [ 1 1 ; 1 2]
%               - Lists of the indices of the selected coeff of the homog tensor *   
%                 I*m : for * written as a matrix
%                 I*v : for * written as a vector
%
% HISTORY
% V. Calisti  01/2021: code implementation.
%**************************************************************************


function [FunJ, GradFunJ, Idx] = init_func(s)   


% Create vC,..,vFh symbolic vector variables of coeff of Ch,...,Fh
% coming in the cost function.
% nC,...,nF are the length of vC,...,vF.
N  = split(s,' ');

%cC = unique(N(find(count(N,'Ch')))) ;  vC = str2sym( cC ) ; nC = length(vC) ;
%cS = unique(N(find(count(N,'Sh')))) ;  vS = str2sym( cS ) ; nS = length(vS) ;
%cD = unique(N(find(count(N,'Dh')))) ;  vD = str2sym( cD ) ; nD = length(vD) ;

cC = unique(N(find(count(N,'Ch')))) ;  vC = sym( cC ) ; nC = length(vC) ;
cS = unique(N(find(count(N,'Sh')))) ;  vS = sym( cS ) ; nS = length(vS) ;
cD = unique(N(find(count(N,'Dh')))) ;  vD = sym( cD ) ; nD = length(vD) ;

% Create arrays lC,...,lFh 
L  = cell2mat( [ cC ; cS ; cD ] );
iJ = [ str2num(L(:,3))  str2num(L(:,4)) ];
lC = []; lS = []; lD = []; 
for i = 1 : nC
    lC(i,:) = iJ(i,:); 
end
for i = 1 : nS
    lS(i,:) = iJ(nC+i,:); 
end
for i = 1 : nD
    lD(i,:) = iJ(nC+nS+i,:); 
end

Idx.lC  = lC ; Idx.lS  = lS ; Idx.lD  = lD ; 

% Create the cost function and its gradient : FunJ and GradFunJ
sVar  = [vC ; vS ; vD];
%f     = str2sym(s);
f     = sym(s); 
Gradf = gradient(f, sVar);

FunJ     = matlabFunction(f, 'Vars', {sVar});
GradFunJ = matlabFunction(Gradf , 'Vars', {sVar});

% load arrays to link indices of tensors following conventions 
% written in docstrings
[I4vm, I5vm, I6vm] = relindice();

% Create indices I*m and I*v for the tensors *
% Chom
if isempty(lC)
    IC = lC ;
else
    for k = 1:size(lC,1)
        iC = find(ismember(I4vm(:,3:4), lC(k,:),'rows'));
        IC(k) = I4vm(iC,1);
    end
    IC = IC';
end
% Shom
if isempty(lS)
    IS = lS ;
else
    for k = 1:size(lS,1)
        iS = find(ismember(I4vm(:,3:4), lS(k,:),'rows'));
        IS(k) = I4vm(iS,1);
    end
    IS = IS';
end
% Dhom
if isempty(lD)
    IDv = lD; IDm = lD; 
else
    for k = 1:size(lD,1)
        iD = find(ismember(I6vm(:,3:4), lD(k,:),'rows'));
        IDm(k) = I6vm(iD,1);
        IDv(k) = I6vm(iD,2);
    end
    IDm = IDm'; IDv = IDv'; 
end

Idx.IC   = IC ;
Idx.IS   = IS ;
Idx.IDm  = IDm;  Idx.IDv  = IDv;

end

