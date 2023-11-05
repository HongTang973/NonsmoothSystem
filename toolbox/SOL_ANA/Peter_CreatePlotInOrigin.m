

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Peter_CreatePlotInOrigin(mdata,data_name)   
    % This will connect to an existing instance of Origin, or create a new one if none exist
    originObj=actxserver('Origin.ApplicationSI');

    % Make the Origin session visible
    invoke(originObj, 'Execute', 'doc -mc 1;');
       
    % Clear "dirty" flag in Origin to suppress prompt for saving current project
    % You may want to replace this with code to handling of saving current project
    invoke(originObj, 'IsModified', 'false');
    
    % Load the custom project CreateOriginPlot.OPJ found in Samples area of Origin installation
    invoke(originObj, 'Execute', 'syspath$=system.path.program$;');
    strPath='';
    strPath = invoke(originObj, 'LTStr', 'syspath$');
    invoke(originObj, 'Load', strcat(strPath, 'Samples\COM Server and Client\Matlab\CreatePlotInOrigin.OPJ'));

    % Create some data to send over to Origin - create three vectors
    % This can be replaced with real data such as result of computation in MATLAB
%     mdata = [0.1:0.1:3; 10 * sin(0.1:0.1:3); 20 * cos(0.1:0.1:3)];
    % Transpose the data so that it can be placed in Origin worksheet columns
%     mdata = mdata';
    % Send this data over to the Data1 worksheet
    invoke(originObj, 'PutWorksheet', 'Data1', mdata);
    
    % Rescale the two layers in the graph and copy graph to clipboard
    invoke(originObj, 'Execute', 'page.active = 1; layer - a');
    invoke(originObj, 'CopyPage', 'Graph1');