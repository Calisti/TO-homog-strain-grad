%**************************************************************************
% SAVE AND STORE RESULTS
%**************************************************************************
%
% HISTORY
% V. Calisti  03/2021: code implementation.
%**************************************************************************


function postproc(cost,matprop,Ch,Dh,infoFin,params,mesh,optionInitsave)   

% ----OPTIONS (Saving options) 

% --Create a latex table with informations
optionTex = 'oui';
%optionTex = 'non';
% --For writing file names begining with the date
optiondate = 'oui';
%optiondate = 'non';
% --For saving init_file.m  'oui'
optioninit = optionInitsave;

% Load the tensors
tensors.Ch = Ch;
tensors.Dh = Dh;

Chstr = matrix2stringtex(tensVoigt(Ch.tens));
Dhstr = matrix2stringtex(tensVoigt(Dh.tens));

% Create a folder with respect to the date : 'the_date'
the_date = datestr(date,'yyyy-mm-dd');
path     = strcat('tests/',the_date);
if ~exist(path, 'dir')
    mkdir(path)
end

% Go to the tests/the_date folder
cd(path)

% Create a sub folder with respect to
% the cost functional: 'CostName'
CostName = num2str(cost.name);
CostName = strrep(CostName,' ', '');
CostName = strrep(CostName,'/', 'on');
CostName = strrep(CostName,'+', 'p');
CostName = strrep(CostName,'*', 'x');
CostName = strrep(CostName,'^', 'E');
CostName = strrep(CostName,'.', '');
CostName = strrep(CostName,'(', '');
CostName = strrep(CostName,')', '');
if ~exist(CostName, 'dir')
    mkdir(CostName)
    % Go to the tests/the_date/Costname folder
    cd(CostName)
    num_files = 0; save('num_files.mat','num_files')
else
    % Go to the tests/the_date/Costname folder
    cd(CostName)
    num_files = load('num_files.mat');
    num_files = num_files.num_files + 1;
    save('num_files.mat','num_files')
end


% For writing files name begining with the date
if strcmp(optiondate, 'oui')
    shortdate = erase(the_date, '-');
    shortdate = shortdate(3:end);
elseif strcmp(optiondate, 'non')
    shortdate = '';
end

% Save the tensors in a .mat file 
str_tens = strcat(shortdate,'homogenized_tensors',num2str(num_files),'.mat');
save(str_tens, 'tensors')

% Save all info about the parameter in a .mat file 
str_info = strcat(shortdate,'info_fin',num2str(num_files),'.mat');
save(str_info, 'infoFin')

% Save the figures : 
%   - Init shape
%   - Final shape
%   - evolution of the shape functional
%   - evolution of theta angle
%   - evolution of the volume
%CostName = strcat(CostName,'r',num2str(ComptRemesh));
Num = num2str(num_files);
saveas(figure(2), strcat(shortdate,'shape_init'  ,Num,'.png') )
saveas(figure(3), strcat(shortdate,'shape_final' ,Num,'.png') )
saveas(figure(4), strcat(shortdate,'shape_func'  ,Num,'.png') )
saveas(figure(5), strcat(shortdate,'angle'       ,Num,'.png') )
% saveas(figure(6), strcat(shortdate,'volume'      ,Num,'.png') )

% Create a latex table with informations
if strcmp(optionTex , 'oui')

    fileTexName = strcat(shortdate,'table',num2str(num_files),'.tex');
    fileID = fopen(fileTexName, 'a');
    % Import the figure
    figName = strcat(shortdate,'shape_final' ,Num,'.png'); %name of fig file
    figTitle= strcat('Final shape for the cost function J=', cost.name);  
    fprintf(fileID, strcat(                                             ...
        '\\begin{minipage}{0.65\\linewidth} \n'                        ,...
        '\\begin{figure}[H] \n'                                        ,...
        '\\fbox{\\includegraphics[width=10cm]{images/',figName,'}} \n'  ,...
        '\\caption{',figTitle,'} \n'                                   ,...
        '\\end{figure} \n'                                             ,...
        '\\end{minipage} \n'                                           ,...
        '\\hfill \n'                                                   ,...
        '\\begin{minipage}{0.25\\linewidth} \n' ) );
    % Write a table containing information about init
    a1 = params.penalty ;
    str1 = strcat('penalty & $' , num2str(a1), '$ \\\\ \n');
    a2 = matprop.gamma ;
    str2 = strcat('$\\gamma$ & $' , num2str(a2), '$ \\\\ \n');
    a3 = mesh.ni ;
    str3 = strcat('initmesh & $' , num2str(a3), '$ \\\\ \n');
    a4 = params.kstart ;
    str4 = strcat('$k_i$ & $' , num2str(a4), '$ \\\\ \n');
    
    fprintf(fileID, '\\begin{tabular}{|c|c|}\n');
    fprintf(fileID, '\\hline   \n');
    fprintf(fileID, str1);
    fprintf(fileID, str2);
    fprintf(fileID, str3);
    fprintf(fileID, str4);
    fprintf(fileID, '\\hline   \n');
    fprintf(fileID, '\\end{tabular} \n');
    fprintf(fileID, '\\end{minipage} \n');
    
    % Write the Chom and Dhom tensors
    fprintf(fileID, '\n');
    fprintf(fileID, '\\begin{align}\n');
    fprintf(fileID, 'C^{hom} & = \n');
    fprintf(fileID, strcat(Chstr, ', \\quad \\theta = ', ...
                    num2str(infoFin.theta(2)) , ' \\\\ \n') );
    fprintf(fileID, 'D^{hom} & = \n');
    fprintf(fileID, strcat(Dhstr, '   \n') );
    fprintf(fileID, '\\end{align}\n');
    
    fclose(fileID);
end


% Save the init file 
if strcmp(optioninit, 'oui')
    copyfile ../../../initialisation/init_file.m
    movefile('init_file.m', strcat('init_file',Num,'.m') )
end

% Go back to dhomo folder
cd ..
cd ../..

end

