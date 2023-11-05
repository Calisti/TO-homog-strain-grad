%**************************************************************************
% CONVERT MATLAB MATRIX INTO TEX ARRAY
%
% Input:  matlab matrix
% Output: string for compiling to tex array
%
% HISTORY
% A.A. Novotny, V. Calisti  2021: code implementation.
%**************************************************************************

function S = matrix2stringtex(A)

leFormat = '% .5G';

S = "\\left( {\\scriptsize \n";
S = strcat(S, '\\begin{array}{');
[I,J] = size(A);
for j=1:J 
    S = strcat(S, 'c');
end
S = strcat(S, "}\n ");

for i=1:I-1
    for j=1:J-1
        s = num2str(A(i,j),leFormat); 
        S = strcat(S, s, " & "); 
    end
    s = num2str(A(i,J),leFormat);
    S = strcat(S, s, " \\\\ \n");  
end

for j=1:J-1
    s = num2str(A(I,j),leFormat);
    S = strcat(S, s, " & ");
end

s = num2str(A(I,J),leFormat);
S = strcat(S, s, " \n");
S = strcat(S, "\\end{array} } \n");
S = strcat(S, '\\right)');

end