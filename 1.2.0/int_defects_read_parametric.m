function [ outputArray ] = int_defects_read_parametric(paramLogStruct, filenamesStruct, varargin)
%int_defects_read_parametric.m
%Harry Coules 2015
%
%DESCRIPTION
%This function is used for reading contour integral results from a large
%number of .dat files created by a parametric series of Abaqus models. It
%is designed to be used as part of the int_defects package.
%
%INPUT ARGUMENTS
%   paramLogStruct - Structure containing information about the parameters
%     used in each finite element model.
%   filenamesStruct - Structure containing filenames.
%   *OPTIONAL*
%   readContoursStr - String defining what contour integral the function
%     should attempt to read. Can contain 'k', 'j' and 't'.
%OUTPUT ARGUMENTS
%   outputArray - An array of structs containing the output data. The heirarchy
%       is: {run no}.integral type{variable no, crack no}(contour no, location no along crack)
%
%% Optional input argument
%This determines what contour integrals the function will attempt to read.
readKFlag = false;
readJFlag = false;
readTFlag = false;
if isempty(varargin)
    %Default behaviour: read K only
    readKFlag = true;
elseif length(varargin) == 1
    if any(strfind(varargin{1},'k')) || any(strfind(varargin{1},'K'))
        readKFlag = true;
    end
    if any(strfind(varargin{1},'j')) || any(strfind(varargin{1},'J'))
        readJFlag = true;
    end
    if any(strfind(varargin{1},'t')) || any(strfind(varargin{1},'T'))
        readTFlag = true;
    end
else
    error('Too many input arguments.');
end

if ~readKFlag && ~readJFlag && ~readTFlag
    error('No contour integral types have been specified to be read!')
end

%% Get directory listing
cd(filenamesStruct.mainFolderName);     %cd to main folder
d = dir;
foldersToIgnore = strcmpi({d(:).name},'.')|...
    strcmpi({d(:).name},'..')|...
    strcmpi({d(:).name},'workspace_dump.mat')|...
    strcmpi({d(:).name},'log1.txt')|...
    strcmpi({d(:).name},'log2.txt')|...
    strcmpi({d(:).name},'matlab_cmd_log1.txt')|...
    strcmpi({d(:).name},'matlab_cmd_log2.txt');
folderNames = {d(~foldersToIgnore).name}';

if length(folderNames)~=length(paramLogStruct)
    warning('The number of folders in main directory is different to the number runs recorded in paramLogStruct.');
end

%% Read .dat files
totalTime = tic;
for k1 = 1:length(folderNames)
    modelTime = tic;
    disp(['READING DATA FROM PARAMETRIC MODEL NO. ',num2str(k1),' OUT OF ',num2str(length(folderNames))]);
    datReadSuccessFlag = false;
    
    %cd to folder
    folderName = char(folderNames(k1));
    cd(folderName);
    
    %Read .dat file
    datFiles = dir('*.dat');
    if isempty(datFiles)
        warning('No valid .dat file detected for this parameter set. Skipping...');
        kOutputArray = {};
        jOutputArray = {};
        tOutputArray = {};
        stepIncTimes = {};
        datReadSuccessFlag = false;
    elseif length(datFiles) == 1
        [~,jobName,~] = fileparts(datFiles.name);
        datScanInpStruct.filename = strcat(jobName,'.dat');
        if readKFlag
            [ kOutputArray, stepIncTimes ] = read_dat_contour_integral( datScanInpStruct , 'k', 1 );  %Note it is assumed that the output for Step 1 is what's required
        end
        if readJFlag
            [ jOutputArray, stepIncTimes ] = read_dat_contour_integral( datScanInpStruct,'j', 1 );
        end
        if readTFlag
            [ tOutputArray, stepIncTimes ] = read_dat_contour_integral( datScanInpStruct , 't', 1 );
        end
        datReadSuccessFlag = true;
    else
        error('There is more than one .dat file in the current Abaqus working folder.');
    end
    
    %Add cell arrays to output array (an array of structs)
    if readKFlag
        outputArray{k1}.k = kOutputArray;
    end
    if readJFlag
        outputArray{k1}.j = jOutputArray;
    end
    if readTFlag
        outputArray{k1}.t = tOutputArray;
    end
    outputArray{k1}.stepIncTimes = stepIncTimes;
    outputArray{k1}.datReadSuccessFlag = datReadSuccessFlag;
    
    %cd back to main directory
    cd ..
    
    %Display
    disp(['Walltime for this model: ',num2str(toc(modelTime),'%4.2f'),' seconds']);
    disp(['Total walltime: ',num2str(toc(totalTime),'%4.2f'),' seconds']);
    disp(' ');
end
cd ..
    
end