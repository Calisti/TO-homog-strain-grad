%**************************************************************************
% Translate a field for having its mean equal to zero  
%**************************************************************************
% DESCRIPTION
% Substract a field with its the avarage, so that the resulting field
% has a null mean.  
%
% INPUT
%  - Ui:      arbitrary field 
%  - psi:     level-set function
%  - mesh:    pdetool mesh struct
%  - matprop: material properties struct
%  - bc:      extra boundry conditions struct
%
% OUTPUT
%  - U:       output field with null mean  
%
% HISTORY
% A.A. Novotny, V. Calisti  05/2020: code implementation.
%**************************************************************************

function U = fixmean(Ui,psi,mesh,matprop,bc)

option = 'nodes' ; 
% option = 'triangle' ; %mauvaise idée
p = mesh.p; t = mesh.t;
np = size(p,2); nU = size(Ui,2) ;
xi = 1 ; xf = np ; yi = np+1 ; yf = 2*np ; 
volRVE = bc.volRVE;   
[~,~,unitF] = assema(p,t,0,0,1);
U = zeros(2*np,nU);
%option = 'do nothing'; U = Ui;

if strcmp(option,'nodes')

    for i=1:nU
        U(xi:xf,i) = Ui(xi:xf,i) - (1/volRVE)*dot(Ui(xi:xf,i),unitF);
        U(yi:yf,i) = Ui(yi:yf,i) - (1/volRVE)*dot(Ui(yi:yf,i),unitF);
    end

end

if strcmp(option,'triangle')
   
    AREA = mesh.area;
    for i=1:nU
        tUx = pdeintrp(p,t,Ui(xi:xf,i));
        tUy = pdeintrp(p,t,Ui(yi:yf,i));
        tUx = tUx - (1/volRVE)*sum(tUx.*AREA);
        tUy = tUy - (1/volRVE)*sum(tUy.*AREA);
        U(xi:xf,i) = pdeprtni(p,t,tUx);
        U(yi:yf,i) = pdeprtni(p,t,tUy);
    end
    
end


end
