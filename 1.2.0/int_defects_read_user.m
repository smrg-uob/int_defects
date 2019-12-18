%int_defects_read_user.m
%Harry Coules 2017
%
%DESCRIPTION
%This is a user-friendly 'wrapper' script for int_defects_read_parametric.
%It automates the file read/write/rename operations that are required for
%using int_defects_read_parametric when there multiple results folders in
%the working directory.
%
%BEHAVIOUR
%int_defects_read_user queries the user for a directory name. The directory
%must contain model results in the standard format used by the int_defects
%package. This script scans the model results files and saves the results
%in *my_abaqus_param_files_directory*\workspace_dump.mat.
%
%% Get directory name from user
clear

%Query user for the name of the directory containing the data
folderName = input('Directory to be read from: abaqus_param_files_','s');
folderName = strcat('abaqus_param_files_',folderName);

%Query the user for the type of contour integral output
contourType = input('Type of contour integral results (k, j or t):   ','s');

%% Read data
movefile(folderName,'abaqus_param_files');  %Rename directory

%Load workspace_dump.mat
cd('abaqus_param_files');
load('workspace_dump.mat');
cd ..

%Read data from .dat files
[ outputArray ] = int_defects_read_parametric(paramLogStruct, filenamesStruct, contourType);

%Save data
cd('abaqus_param_files');
save('workspace_dump.mat');
cd ..

movefile('abaqus_param_files',folderName);  %Rename abaqus_param_files_... back to its original name.

disp('End of int_defects_read_user.')