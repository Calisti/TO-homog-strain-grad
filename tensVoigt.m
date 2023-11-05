%**************************************************************************
% Computes the Voigt homogenized tensor
%**************************************************************************
% DESCRIPTION
% Rewrite the input homogenized T tensor of order 4 or 6
% with the Voigt convention --> [T]
% A   = A_ij (ei x ej)        --->
% [A] = [ A11 A22 sqrt(2)A12 ]'
% B   = B_ijk (ei x ej x ek)  ---> 
% [B] = [ B111 B221 sqrt(2)B121 B112 B222 sqrt(2)B122 ]'
% Thus :
% A:T:A = [A]:[T]:[A]  and  B:T:B = [B]:[T]:[B]
%
% INPUT
%  - T:       input homogenized tensor
%
% OUTPUT
%  - Tv:      homogenized tensor written in Voigt convention : Tv = [T]
%
% HISTORY
% V. Calisti  02/2021: code implementation.
%**************************************************************************

function Tv = tensVoigt(T)

i = size(T);

MatCofVoigt = [   1        1      sqrt(2) ;...
                  1        1      sqrt(2) ;...
                sqrt(2)  sqrt(2)     2      ];

if i(1) == 3
    
    Tv = MatCofVoigt .* T ;
    
elseif i(1) == 6
    
    % Base Voigt vectors column written in init base  
    MatBaseChg = [ 1 0 0 0 0 0 ;...
                   0 0 0 1 0 0 ;...
                   0 1 0 0 0 0 ;...
                   0 0 0 0 1 0 ;...
                   0 0 1 0 0 0 ;...
                   0 0 0 0 0 1     ];

    Tv = repmat(MatCofVoigt, 2, 2) .* ( (MatBaseChg') * T * MatBaseChg );
    
else
    
    error('Size Error: the input tensor of tensVoigt method is not of order 4 or 6')
   
end

end