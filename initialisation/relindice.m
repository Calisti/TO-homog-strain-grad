%**************************************************************************
% Indices relation
%**************************************************************************
% DESCRIPTION
% Compute matrices of indices used for the identification of the 
% coefficients of the homogenized tensors T and of their topological
% derivatives DtT. These matrices are called I*vm, and are composed with 
% four columns, and as much rows as the number of the coefficients of 
% tensors of order *. These matrices are to be read row by row. 
% For a given tensor T, for a given row [a b i j], we have:  
% -- [i j] are the indices of the coeff (T_ij)_ij when T is written in
%    our matrix convention 
% -- [a] is the index such that T(a) corresponds to T_ij
% -- [b] is the index such taht DtT(b) corresoinds to DtT_ij
%
% OUTPUT
%  - I4vm : index matrice for 4th order tensors
%  - I5vm : index matrice for 5th order tensors
%  - I6vm : index matrice for 6th order tensors
%
% HISTORY
% V. Calisti  01/2021: code implementation.
%**************************************************************************

function [I4vm, I5vm, I6vm] = relindice()


I4vm = [1 1 1 1 ; ...
        4 2 1 2 ; ...
        7 3 1 3 ; ... 
        5 4 2 2 ; ...
        8 5 2 3 ; ...
        9 6 3 3      ];
    
I5vm = [1   1   1 1 ; ...
        4   2   1 2 ; ...
        7   3   1 3 ; ... 
        10  4   1 4 ; ...
        13  5   1 5 ; ...
        16  6   1 6 ; ... 
        2   7   2 1 ; ...
        5   8   2 2 ; ...
        8   9   2 3 ; ... 
        11  10  2 4 ; ...
        14  11  2 5 ; ...
        17  12  2 6 ; ... 
        3   13  3 1 ; ...
        6   14  3 2 ; ...
        9   15  3 3 ; ... 
        12  16  3 4 ; ...
        15  17  3 5 ; ...
        18  18  3 6 ];

    
I6vec = (1:36)';
I6mat = zeros(36,2);
A = zeros(6,6); A(:) = 1:36 ; I = zeros(1,36);

t1 = 1;
for i = 1:3
    for j = 1:3
        %
        I6mat(t1,:) = [2*i-1,2*j-1] ;
        I(t1) = A(2*i-1,2*j-1) ;
        t1 = t1 + 1;
        %
        I6mat(t1,:) = [2*i  ,2*j  ] ;
        I(t1) = A(2*i  ,2*j  ) ;
        t1 = t1 + 1;
        %
        I6mat(t1,:) = [2*i-1,2*j  ] ;
        I(t1) = A(2*i-1,2*j  ) ;
        t1 = t1 + 1;
        %
        I6mat(t1,:) = [2*i  ,2*j-1]  ;
        I(t1) = A(2*i  ,2*j-1) ;
        t1 = t1 + 1; 
    end
end

I6vm = [I(1:36)' , I6vec , I6mat];

end