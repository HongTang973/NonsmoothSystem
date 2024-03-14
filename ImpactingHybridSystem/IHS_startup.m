%% -------------  detect the platform where the code is runing -------- %%
current_folder = pwd;
isUoB = get_platform();
%>
if isUoB; onedrive_path = 'C:\Users\ib20968\'; else; onedrive_path = 'F:\onedrive\'; end
addpath([onedrive_path, 'OneDrive - University of Bristol\Codes stall\sys'])
%> include the general toolbox for customized commands
addpath(genpath([onedrive_path,'OneDrive - University of Bristol\Codes stall\General_Toolbox']))

%> include the continuation tool box
cd([ onedrive_path,'OneDrive - University of Bristol\Codes stall\coco_2020Mar22\coco'])
startup;
%> include the path of the toolbox for the IHS 
addpath(genpath([onedrive_path,'OneDrive - University of Bristol\Codes stall\NonsmoothSystem\']))
%>
cd(current_folder)

