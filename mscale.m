%**************************************************************************
% Micro-Structure Topology Design Optimization with Volume Constraint 
%**************************************************************************
%
% DESCRIPTION
% Computes the topological derivative and use it together with a level-set 
% domain representation method in the topology design of micro-structutes
%
% HISTORY
% S. Amstutz & A.A. Novotny      06/2009: code implementation.
% S.M. Giusti                    03/2011: code updating.
% A.A. Novotny                   06/2012: code updating.
% V.   Calisti                   2020-21: code updating.
%**************************************************************************

%%%% FOR WINDOWS
%clear all; close all; dir = 'initialisation'; 
% load problem data
%cd('initialisation')
%     [mesh, pdecoef, matprop, params, bc, psi0, cost] = init_file;
%cd ..

%%%% FOR LINUX
clear; close all; dir = 'initialisation';    
% Load problem data
cd('initialisation')
    [mesh, pdecoef, matprop, params, bc, psi0, cost] = init_file; 
cd ..

% Mesh and geometry  parameter
p = mesh.p; e = mesh.e; t = mesh.t; g = mesh.g; area = mesh.area;
count_remesh = 0; 
% Set multi-scale dof's and boundary condition parameters
[mesh,bc] = splitmesh(mesh,size(pdecoef.f,1),bc); volRVE = bc.volRVE;
% Topology optimization parameters
stop = params.stop; kmin = params.kmin; penalization = params.penalization; 
volfrac = params.volfrac; auglag = params.auglag; voleps = params.voleps;
% Freeze nodes 
phold = [];
for i = 1:size(mesh.ghold,1)
    phold = cat(1,phold,unique(t(1:3,t(4,:)==mesh.ghold(i))));
end
% free nodes
pfree = setdiff(1:size(p,2),phold);
% mass matrix of unity density
[~,unitM,~] = assema(p,t,0,1,0);       
% level-set function nomalization 
psi = psi0/sqrt(dot(unitM*psi0,psi0)); 

% Compute the correctors and homogenized tensors 
[U,Ut,K] = solve(psi, mesh, matprop, pdecoef, bc);
Ch0 = tens_Ch(mesh,U,matprop,psi,volRVE);
Utt = solve2(U,Ut,Ch0,K,psi, mesh, matprop, bc,volRVE);
Dh0 = tens_Dh(U,Ut,Utt,Ch0,mesh,matprop,psi,volRVE);
% Save the initial tensors in Th.init (used for Func?)
tens_init.Ch = Ch0;
tens_init.Dh = Dh0;
Th.tinit     = tens_init;                     

% Shape functional associated to the initial guess
params.sfJ0 = abs( shfuncJ(cost,Ch0,Dh0) );
%=================================
%params.sfJ0 = 1
%=================================
% volume associated to the initial guess
% characteristic function
tchi = pdeintrp(p,t,(psi0 < 0)); 
% update the volume of the bulk phase
volinit = dot(area,tchi),       
%=================================
%volinit = 1
%=================================
% initial volume
params.volinit = volinit; 
voltarget = volinit * volfrac;

% Stop criterion over the volume constraint 
volstop = voleps * volinit;
if penalization == 1
    volstop = volinit;
end

% Plot hold-all domain
figure(2); clf;                                                            
set(2,'WindowStyle','docked');                                          
pdeplot(p,e,t,'xydata',(psi>0),'xystyle','flat','colormap','gray',...
              'xygrid','off','colorbar','off'); axis image; axis off;% pause
            
          
% Initializing time parameter 'k' of the line search 
kstart = params.kstart; 
k = kstart; 
iter = 0; 
% Initializing lists of results
% Value of sh. func (gsf), angle (gth) and volume (gvol) 
gsf = []; gth = []; gvol = []; 
% Save for iteration (I) when change of constrast (GAM), remesh (M),
% shape functional (S), kstart (KS)
I=[0;0]; GAM=[0;0]; M=['.';'.']; S={'.';'.'}; KS = [0;0];

option = 'null';

while not(strcmp(option,'s'))
    
    iter = iter + 1;
    if(iter == 200) 
        option = 's';
    end
    
    % Characteristic function.
    tchi = pdeintrp(p,t,(psi < 0)); 
    % Update the volume of the bulk fase.
    volume = dot(area,tchi);
    % Compute correctors and homogenized tensors.
    [U,Ut,K] = solve(psi, mesh, matprop, pdecoef, bc);
    Ch = tens_Ch(mesh,U,matprop,psi,volRVE);
    Utt = solve2(U,Ut,Ch,K,psi, mesh, matprop, bc,volRVE);
    Dh = tens_Dh(U,Ut,Utt,Ch,mesh,matprop,psi,volRVE);
    % Compute Shape Functional.
    sf = shfunc(cost,Ch,Dh,volume,params);
    % Compute function g: g = -DT (bulk) and g = DT (inclusion)   
    dt = topder(mesh,U,Ut,Utt,K,bc,volume,matprop,psi,params,volRVE,cost,Ch,Dh);
    dt = dt/sqrt(dot(unitM*dt,dt));  % g function normalization
    dt(phold) = psi(phold);          % freeze dt function
    cosin = max(min(dot(unitM*dt,psi),1.0),-1.0);
    theta = max(real(acos(cosin)),1.0e-4);
       
    % Performs a line-search 
    sfold = sf; psiold = psi; sf = sf + 1; k = min(kstart,2*k); 
    while  (sf > sfold ) && (k > kmin/4)
        % Update level-set & charct functions and volume of the bulk fase
        psi(pfree)  = (sin((1-k)*theta)*psiold(pfree)...
                    +  sin(k*theta)*dt(pfree))./sin(theta);
        %psi = psi/sqrt(dot(unitM*psi,psi)); %normalization
        tchi = pdeintrp(p,t,(psi < 0));
        volume = dot(area,tchi);     
        % Update correctors and homog. tensors
        [U,Ut,K] = solve(psi,mesh,matprop,pdecoef,bc); 
        Ch = tens_Ch(mesh,U,matprop,psi,volRVE);
        Utt = solve2(U,Ut,Ch,K,psi, mesh, matprop, bc,volRVE);
        Dh = tens_Dh(U,Ut,Utt,Ch,mesh,matprop,psi,volRVE);
        % Update shape func.
        sf = shfunc(cost,Ch,Dh,volume,params);
        k = k / 2;               
    end   
    psi = psi/sqrt(dot(unitM*psi,psi)); %normalization
    
    k = k * 2;
    
    gsf  = [gsf,sf];
    gth  = [gth,theta];
    gvol = [gvol,volume*100/volinit];

    disp('Optimization procedure'); 
    disp('----------------------'); 
    disp(['iter        = ', num2str(iter)]); 
    disp(['volume      = ', num2str(volume),' => ', ...
                            num2str(volume*100/volinit), '%']);
    disp(['sf          = ', num2str(sf)]);     
    disp(['k           = ', num2str(k)]);          
    disp(['theta       = ', num2str(theta*180/pi)]);  
    disp(['penalty     = ', num2str(params.penalty)]);     
   
    drawnow  
    % Black and white shape
    figure(3); clf;
    set(3,'WindowStyle','docked');  
    pdeplot(p,[],t,'xydata',(psi>0),'xystyle','flat','colormap','gray',...
                  'xygrid','off','colorbar','off'); axis image; axis on;
    %hold on ; pdeplot(p,e,t); hold off % +MESH
    % Shape func evolution
    figure(4); clf; plot(gsf, ':*k'); 
    title('Shape Function'); 
    set(4,'WindowStyle','docked');   
    % Theta evolution
    figure(5); clf; plot(gth*180/pi, ':*k'); 
    title('Theta Angle'); 
    set(5,'WindowStyle','docked'); 
    % Volume evolution
    figure(6); clf; plot(gvol, ':*k'); 
    title('Volume'); 
    set(6,'WindowStyle','docked'); 
    % Level-set : Psi
    figure(7);clf;
    pdeplot(p,e,t,'xydata',psi,'Contour','on'); axis image; %axis off;
    set(7,'WindowStyle','docked'); 
    %hold on; pdeplot(p,e,t); hold off; % +MESH
    % Signed topological derivative : g
    figure(8);clf;
    pdeplot(p,e,t,'xydata',dt,'Contour','on'); axis image; %axis off;
    set(8,'WindowStyle','docked'); 
    %hold on; pdeplot(p,e,t); hold off; % +MESH
    
    if or(k < kmin, theta < stop)
        
        % The volume constraint is not satisfies (for penalization==2 or 3)
        difvol = volume-voltarget;
        if (abs(difvol) > volstop)  %&& (k > 0.00001)%######
            disp('difvol')
            if penalization == 2 % increase the penalization parameter
               params.penalty = 2.0 * params.penalty;
               k = kstart;
            elseif penalization == 3 % increase the lagrangian multiplier
                coef = volume / voltarget; tau = auglag;
                penalty = params.penalty;
                penalty = penalty + tau/auglag ...
                        * (max(0,penalty - auglag*(1.0-coef)) - penalty);
                params.penalty = penalty; 
                k = kstart;
            end
            
        % Stop or try modifications 
        elseif or( (abs(volume-voltarget) <= volstop) , k < kmin )            
           
            option = 'null';                                              
            while not(strcmp(option,'m')) && not(strcmp(option,'s')) ...
                                          && not(strcmp(option,'c'))  
                str = '\n -> Enter "m" to modify or "s" to stop:  \ ';
                option = input(str, 's');
            end   
            
            % Modifications:
            if (option == 'm')
                
                optionR = 'null'; optionP = 'null'; 
                I = [I; iter]; 
                % Choice: homogeneous remesh, local remesh, or no remesh 
                while not(strcmp(optionR,'h')) && not(strcmp(optionR,'j'))...
                                               && not(strcmp(optionR,'n')) 
                    str = '\n -> Enter "h" for homog. remesh, "j" for local remesh, or "n" for no remesh: \ ';
                    optionR = input(str, 's'); 
                    M = [M; optionR];
                end
                % Choice: change some parameters or not
                while not(strcmp(optionP,'p')) && not(strcmp(optionP,'n')) 
                    str = '\n -> Enter "p" to change params, or "n" for no change: \ ';
                    optionP = input(str, 's');
                end
                
                % Remesh:
                %--------
                if strcmp(optionR,'j')
                    % Local remesh
                    count_remesh = count_remesh + 1;% Count numb. of remesh 
                    [p,e,t,psi,psi0] = remeshlocal(mesh,psi,psi0);
                    mesh.p = p; mesh.e = e; mesh.t = t;
                elseif strcmp(optionR,'h')                    
                    % Homogeneous remesh
                    count_remesh = count_remesh + 1;% Count numb. of remesh 
                    for i = 1:(strcmp(mesh.remesh,'longest')*2 + strcmp(mesh.remesh,'regular'))
                        [p,e,t,psi] = refinemesh(g,p,e,t,[psi psi0],mesh.remesh);
                        % psi is the updated level set, and psi0 is used
                        % for updating the initial volume
                        psi0 = psi(:,2); psi = psi(:,1);
                        mesh.p = p; mesh.e = e; mesh.t = t;
                    end
                end
                % Update mesh parameters
                [~,unitM,~] = assema(p,t,0,1,0);
                area = pdetrg(p,t); mesh.area = area; k = kstart;
                [mesh,bc] = splitmesh(mesh,size(pdecoef.f,1),bc);
                % Freeze groups
                phold = [];
                for i = 1:size(mesh.ghold,1)
                    %nodes to be fixed
                    phold = cat(1,phold,unique(t(1:3,t(4,:)==mesh.ghold(i)))); 
                end
                pfree = setdiff(1:size(p,2),phold); 
                option = 'null';
                                               
                tchi0 = pdeintrp(p,t,(psi0 < 0)); 
                volinit = dot(area,tchi0); 
                %=================================
                %volinit = 1;
                %=================================                
                params.volinit = volinit; 
                voltarget = volinit * volfrac;
                % stop criterion over the volume constraint 
                volstop = voleps * volinit;
                if penalization == 1
                    volstop = volinit;
                end
                
                % Change parameters:
                %-------------------
                if strcmp(optionP,'p')
                    [cost,matprop,params.sfJ0,params.kstart] = ...
                        modifParams(mesh,psi,matprop,pdecoef,  ...
                                        bc,volRVE,cost,params);  
                    kstart = params.kstart; 
                    k = kstart; 
                end
                GAM = [GAM; matprop.gamma]; S{size(S,1)+1,1} = cost.name;
                KS  = [KS; kstart];
                
            end
            
        end
        
    end
    
end
              
fig = input(' -> Export figures? "y" or "n" : ', 's');           

if strcmp(fig,'y')

    % Load chosen information to save 
    infoFin.iter  = iter;   infoFin.vol  = volume;
    infoFin.sf    = sf;     infoFin.k    = k; 
    infoFin.pen   = params.penalty;
    infoFin.func  = cost.name;
    infoFin.theta = [theta, theta*180/pi]; 
    infoFin.gsf   = gsf; infoFin.gth = gth; infoFin.gvol = gvol;
    infoFin.modif = table(I,M,GAM,S,KS);
    %=================== Save and store the figures
    postproc(cost,matprop,Ch,Dh,infoFin,params,mesh,'oui') 
    %=================== 
    
else
    
    return ;
    
end


