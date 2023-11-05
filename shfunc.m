%**************************************************************************
% Shape Functional Associated to the Compliance with Volume Constraint 
%**************************************************************************
% DESCRIPTION
% Computes the shape functional
% 
% INPUT
%  - cost:   option for the choice of the shape functional
%  - Ch,Dh:  homogenized tensors
%  - volume: current volume
%  - params: topology optimization parameters struct
%
% OUTPUT
%  - sf: shape functional value
%
% HISTORY
% S. Amstutz & A.A. Novotny      06/2009: code implementation.
% J-M.C. Farias                  12/2010: code updating.
% D.E. Campe�o   
% A.A. Novotny
% S.M. Giusti                    03/2011: code updating.
%**************************************************************************

function sf = shfunc(cost,Ch,Dh,volume,params)

    penalty = params.penalty; volinit = params.volinit; auglag= params.auglag;
    volfrac = params.volfrac; voltarget = volinit*volfrac;
    %normalizing shape functional or ntot
    sf = shfuncJ(cost,Ch,Dh) / (1 + cost.normzn*(params.sfJ0 - 1)) ;  
    
    if params.penalization == 1 % linear penalization
        sf = sf + penalty * volume / volinit;
    elseif params.penalization == 2 % exact quadratic penalization
        coef = volume / voltarget;
        sf = sf + penalty * ((coef <= 1)*0.0+(coef > 1)*(1.0-coef)^2);
    elseif params.penalization == 3 % augmented lagrangian multiplier
        coef = volume / voltarget;
        aux = max(coef-1.0, -penalty / auglag);
        sf = sf + penalty * aux + (aux)^2 * auglag / 2.0; 
    end

end

