%**************************************************************************
% Remesh + Change of parameter 
%**************************************************************************
% DESCRIPTION
% Allows to change the different parameters of the problem during the 
% optimization process.
%
% INPUT
% All information about the geometry and properties of the problem
%
% OUTPUT
% Updated : - criterion (cost) and its new init. value (sfJ0)
%           - material constrast property (matprop)
%           - step paramater (k)
%
% HISTORY
% V. Calisti  04/2021: code implementation.
%***********************************************************************

function [cost,matprop,sfJ0,kstart] = modifParams(mesh,psi,matprop,pdecoef,bc,volRVE,cost,params)

% Change parameters 
% -----------------

% * Contrast gamma
str1 = strcat('\n current contrast: gamma = \ ',num2str(matprop.gamma));
str2 = '\n -> enter new contrast: gamma = \ ';
matprop.gamma = input(strcat(str1,str2)); 

% * time-step kstart
str1 = strcat('\n current time-step: kstart = \ ',num2str(params.kstart)); %NEW
str2 = '\n -> enter new time-step: kstart = \ '; %NEW
kstart = input(strcat(str1,str2)); %NEW

% % * Shape Function
% str1 = strcat('\n current shape function: s = \ ',num2str(cost.name));
% str2 = '\n -> enter new shape function: s = \ ';
% cost.name = input(strcat(str1,str2)); 
% cd('initialisation')
% [FunJ, GradFunJ, Idx] = init_func(cost.name);
% cost.FunJ             = FunJ;
% cost.GradFunJ         = GradFunJ;
% cost.Idx              = Idx;
% cd ../..

% Compute the needed homogenized tensors and correctors
[U,Ut,K] = solve(psi, mesh, matprop, pdecoef, bc);
Ch0 = tens_Ch(mesh,U,matprop,psi,volRVE);
Utt = solve2(U,Ut,Ch0,K,psi, mesh, matprop, bc,volRVE);
Dh0 = tens_Dh(U,Ut,Utt,Ch0,mesh,matprop,psi,volRVE);
Th0 = tens_all(cost,U,Ut,Utt,mesh,matprop,psi,volRVE); %NEW
                        
% New init shape functional evaluation after parameters changes 
sfJ0 = abs( shfuncJ(cost,Ch0,Dh0,Th0) );

end

