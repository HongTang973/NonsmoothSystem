%% -------------  detect the platform where the code is runing -------- %%
current_folder = pwd;
[~, name] = system('hostname');
%> the name is with a /n, delete it; in matlab line feed's ASCII is 10
if name(end)== 10; name(end) = []; else; end
%
if strcmp(name,'IT079952')%C:\Users\ib20968\OneDrive
    isUoB = true;
elseif strcmp(name,'LAPTOP-47T2DBU0') %> personal computer
    isUoB = false;
else
    error('Unkonw computer is now owned by you!')
end
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

