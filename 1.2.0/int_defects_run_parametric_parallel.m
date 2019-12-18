function int_defects_run_parametric_parallel( paramLogStruct, filenamesStruct, abqVer, timeOutTime, noWorkers, noCPUs, varargin)
%int_defects_run_parametric_parallel.m
%Harry Coules 2015
%
%DESCRIPTION
%This function is used to run a parametric series of Abaqus models
%representing structures containing crack-like defects. It is used as part
%of the int_defects software package to study defect interaction.
%
%TOOLBOX REQUIREMENT
%This code is designed to use the Parallel Computing Toolbox to run several
%Abaqus jobs in parallel. However, if the PCT is not available then it will
%default to executing one model at a time and the input argument noWorkers
%will be ignored.
%
%INPUT ARGUMENTS
%   filenamesStruct - Structure containing filenames relavent to analysis.
%       .mainFolderName - Name of a main folder containing numbered
%           subfolders which themselves contain the analysis files.
%       .filenameGeneratedInp
%   abqVer - Version of Abaqus/Standard to use. Should normally be given as
%       'abaqus', in which case the default installed Abaqus version will
%       be used. However, 'abq6121' or 'abq6141' can be given to force
%       model execution using older Abaqus versions.
%   timeOutTime - Time (in seconds) after which to kill an Abaqus job if it
%       has not completed. This may be given as Inf if no model time limit
%       is required.
%   noWorkers - Number of parallel MATLAB workers to use. If the Parallel
%       Computing Toolbox is not available, this argument will be ignored.
%   noCPUs - Number of CPUs to instruct Abaqus to use. If this is given as
%       NaN, the number of CPUs is not specified in the command used to run
%       the Abaqus job.
%
%*OPTIONAL*
%   abaqusCleanupFlag - A logical to indicate whether the files that Abaqus
%       generates (except the .dat file) should be deleted after each run 
%       to save space. If not specified, the user will be asked about this
%       when the script is run.
%   contDepCheckFlag - A logical to indicate whether or not
%       contour-dependence checks should be performed after the main run of
%       models. It may be desirable to prevent this, eg. if performing
%       J-integral models which are going to be post-processed manually.
%
% nB. If contDepCheckFlag is given, it must be given in addition to
%   abaqusCleanupFlag. The order of optional input arguments must be:
%   1. abaqusCleanupFlag, 2. contDepCheckFlag.
%
%OUTPUT ARGUMENTS
%   %NONE%
%
%OUTPUT FILES
% - Results files from the Abaqus models are created in sub-directories
%   within the directory abaqus_param_files.
% - workspace_dump.mat is created in the directory abaqus_param_files.
% - Logfiles are created: matlab_cmd_log1.txt, matlab_cmd_log2.txt and (in
%   abaqus_param_files) log1.txt.
%
%
%DEVELOPMENT NOTES
% 12/2017 - Modified to use parpool or matlabpool depending on the
%   available Matlab version.
% 10/2018 - Added the contDepCheckFlag optional input argument and
%   associated functionality.
% 10/2018 - Added read-back of parallel run parameters and function status
%   display.
%
%% Preliminaries
%Get directory listing
cd(filenamesStruct.mainFolderName);
d = dir;
foldersToIgnore = strcmpi({d(:).name},'.')|...
    strcmpi({d(:).name},'..')|...
    strcmpi({d(:).name},'workspace_dump.mat')|...
    strcmpi({d(:).name},'log1.txt');   %Note specific exemptions for files which may also be in there.
folderNames = {d(~foldersToIgnore).name}';
if length(folderNames) ~= length(paramLogStruct);
    warning('Unexpected number of items in working directory.')
    error('Number of folders in the working directory is different from the number of entries in paramLogStruct loaded from workspace_dump.mat.');
end

%Optional input arguments
if length(varargin) > 2
    error('Too many input arguments.');
elseif length(varargin) == 2
    abaqusCleanupFlag = varargin{1};
    contDepCheckFlag = varargin{2};
elseif length(varargin) == 1
    abaqusCleanupFlag = varargin{1};
    contDepCheckFlag = true;
else
    contDepCheckFlag = true;
end

%If abaqusCleanupFlag was not specified, ask the user whether .odb files
%should be deleted.
if ~exist('abaqusCleanupFlag','var')
    abaqusCleanupString = input('Run abaqus_cleanup after each model, so that only .dat and .sta files are saved? (y/n) [default: y]:   ','s');
    if isempty(abaqusCleanupString)
        abaqusCleanupFlag = true;
    elseif strcmpi(abaqusCleanupString,'y')
        abaqusCleanupFlag = true;
    elseif strcmpi(abaqusCleanupString,'n')
        abaqusCleanupFlag = false;
    else
        error('Input should be either ''y'' or ''n''.');
    end
end

%Check whether the parallel computing toolbox is available
parallelAvailFlag = ~isempty(ver('distcomp'));
if ~parallelAvailFlag
    warning('Parallel Computing Toolbox is not available. Defaulting to a single worker. The requested number of parallel jobs will be ignored.');
end

%Read back basic information about the parallel run
disp(' ');
disp('READY TO RUN THE PARAMETRIC SET OF MODELS.')
if ~parallelAvailFlag
    disp(['The number of parallel workers is: ',num2str(noWorkers),'.']);
else
    disp('MATLAB parallelisation is not being used.');
end
disp(['The number of CPUs per worker is: ',num2str(noCPUs),'.']);
disp(['The timeout for each analysis is: ',num2str(timeOutTime,'%.0f'),' sec. (i.e. ',num2str(timeOutTime/60,'%.2f'),' min, or ',num2str(timeOutTime/3600,'%.2f'),' hr)']);
disp(['The version of Abaqus/Standard that will be used by preference is: ',abqVer]);
if abaqusCleanupFlag   %Abaqus output file cleanup
    disp('Abaqus output files (excluding .dat) will be DELETED after each model is run.');
else
    disp('Abaqus output files will be PRESERVED after each model is run.');
end
if contDepCheckFlag   %Contour-dependence checking
    disp('Contour-dependence checks WILL be performed on contour integral output.');
else
    disp('Contour-dependence checks WILL NOT be performed on contour integral output.');
end
disp(' ');
disp('NOTE: Messages from Abaqus/Standard in the following section may appear out-of-order due to parallelisation.');
disp(' ');

%% Run Abaqus analyses
%Create pool of parallel workers
warning('off','MATLAB:datetime:NonstandardSystemTimeZone');
if parallelAvailFlag
    if verLessThan('matlab','8.4')
        matlabpool(noWorkers);
    else
        parPool1 = parpool(noWorkers);
    end
end
if license('test','Distrib_Computing_Toolbox')
    pctRunOnAll warning('off','MATLAB:lang:cannotClearExecutingFunction');
end

%Run using parfor
totalTime = tic;
jobExitStatusSliced = zeros(length(folderNames),1);
parfor k1 = 1:length(folderNames)
    %Preallocate temporaty variables
    jobExitStatus = 0;
    noFolders = 0;
    
    %Run
    modelTime = tic;
    disp(['RUNNING PARAMETRIC MODEL NO. ',num2str(k1),' OUT OF ',num2str(length(folderNames))]);
    
    %cd to folder
    folderName = char(folderNames(k1));
    cd(folderName);
    
    %Check if there is a model in this folder which needs to be run
    runThisModel = true;
    if paramLogStruct(k1).inpStatus
        if logical(length(dir('*.inp'))) && ~logical(length(dir('*.dat')))
            %No action
        elseif ~logical(length(dir('*.inp'))) && logical(length(dir('*.dat')))
            %.dat only
            warning('Model in this folder appears to have run already. Skipping this model.');
            runThisModel = false;
        elseif logical(length(dir('*.inp'))) && logical(length(dir('*.dat')))
            %Both .inp and .dat
            warning(sprintf('Both .inp and .dat files detected in working folder.\nMay be the result of a prior model run.\nAttempting to delete old files and re-run...'));
            abaqus_cleanup({'.com','.dat','.log','.msg','.odb','.prt','.sim','.sta'});
        else
            %Neither .inp nor .dat
            error('Neither .inp nor .dat files exist in this folder')
        end
    else
        warning('No valid .inp file available for this parameter set. Skipping...');
        runThisModel = false;
    end
    
    %Run job (if required)
    if runThisModel
        %Get name of the job
        [~,jobName,~] = fileparts(filenamesStruct.filenameGeneratedInp);
        
        %Run job
        [jobExitStatus] = abaqus_run(jobName,timeOutTime,noCPUs,false,abqVer);
        jobExitStatusSliced(k1) = jobExitStatus; %parfor sliced output variable used for returning the exitStatus of each job
        
        %Clean up larger Abaqus files if required... But not the inp, since this may be
        %needed later if the model fails contour-dependence checks
        if abaqusCleanupFlag
            abaqus_cleanup({'.odb','.prt','.filepart'});
        end
    else
        disp(['Not necessary to run model no. ',num2str(k1),'. runThisModel = 0.']);
        jobExitStatusSliced(k1) = 2;    %Status of 2 indicates that it was not necessary to run this model.
    end
    
    %Return to main directory
    cd ..
    
    %Check how many .dat files have been written so far and how many total
    %folders there are.
    if mod(k1,10)==0    %Only do this for every 10th model, to save time.
        %Note - Since requires a recursive dir, it will take a significant
        %amount of time when there are many model folders.
        dPar = dir;
        foldersToIgnorePar = strcmpi({dPar(:).name},'.')|...
            strcmpi({dPar(:).name},'..')|...
            strcmpi({dPar(:).name},'workspace_dump.mat')|...
            strcmpi({dPar(:).name},'log1.txt')|...
            strcmpi({dPar(:).name},'log2.txt');
        noFolders = length(dPar)-sum(foldersToIgnorePar);
        if ispc
            noDatFiles = length(rdir('*\*.dat'));
        else
            noDatFiles = length(rdir('*/*.dat'));
        end
    else
        noDatFiles = NaN;
    end
    
    %Display
    if runThisModel
        if jobExitStatus == 1
            disp('DONE - a model was run for this set of parameters.');
        else
            disp('FAILED - a job was submitted but failed to run.')
        end
    else
        disp('NO JOB ATTEMPTED - no model was run for this set of parameters.');
    end
    disp(['Walltime for this model: ',num2str(toc(modelTime),'%4.2f'),' seconds']);
    disp(['Total walltime: ',num2str(toc(totalTime),'%4.2f'),' seconds']);
    if ~isnan(noDatFiles)
        disp([num2str(noDatFiles),' .dat files written so far, in ',num2str(noFolders),' folders.']);
    end
    disp(' ');
end
if parallelAvailFlag
    if verLessThan('matlab','8.4')
        matlabpool close;
    else
        delete(parPool1);
    end
end
disp(' ');
disp('FINISHED RUNNING THE PARAMETRIC SET OF MODELS.')
disp(' ');

%% Check for contour-dependent results
%Note that a for loop rather than a parfor loop is used for this. This is
%because to check the result Matlab needs to access the previously-written
%.dat file for each model, and this can cause problems for parfor.
contDepWarnFlagLevel1 = [];
contDepWarnFlagLevel2 = [];
for k1 = 1:length(folderNames)
    if contDepCheckFlag
        %Display
        disp(['Checking contour integral results in folder ',char(folderNames(k1)),' for path independence.']);
        
        %cd to folder
        folderName = char(folderNames(k1));
        cd(folderName);
        
        if logical(length(dir('*.dat')))    %If there's a .dat file in the folder for this model.
            %Get name of the job
            datFiles = dir('*.dat');
            if length(datFiles) > 1
                error('There is more than one .dat file in the current Abaqus working folder.');
            end
            [~,jobName,~] = fileparts(datFiles.name);
            
            %Determine the number of contour integral output requests for this job
            noContourRequests = 0;
            if isfield(paramLogStruct(k1).naturalParams.outputRequests,'contourNameA')
                noContourRequests = noContourRequests + 1;
            end
            if isfield(paramLogStruct(k1).naturalParams.outputRequests,'contourNameB')
                noContourRequests = noContourRequests + 1;
            end    
            
            %Check the contour-dependence. NOTE: It is assumed that for multiple cracks, the type of contour integral output requested is the same.
            if strcmpi(paramLogStruct(k1).naturalParams.outputRequests.contourTypeA,'K_FACTORS')
                [ kCheckArray ] = read_dat_contour_integral( [jobName,'.dat'], 'k', 1);
                [ contDepFlags ] = int_defects_contour_independence_check( kCheckArray, noContourRequests, 'k', 0.01 );
            elseif strcmpi(paramLogStruct(k1).naturalParams.outputRequests.contourTypeA,'J_INTEGRAL')
                [ jCheckArray ] = read_dat_contour_integral( [jobName,'.dat'], 'j', 1);
                [ contDepFlags ] = int_defects_contour_independence_check( jCheckArray, noContourRequests, 'j', 0.03 );    %Note greater tolerance for J then for K
            elseif strcmpi(paramLogStruct(k1).naturalParams.outputRequests.contourTypeA,'T_STRESS')
                warning('Contour-dependence checking not yet implemented for T-stress results.');
                contDepFlags = false;
            elseif strcmpi(paramLogStruct(k1).naturalParams.outputRequests.contourTypeA,'none')
                warning('Contour-dependence checking requested, but no contour integral output has been requested for Crack A. Returning contDepFlags = false.');
                contDepFlags = false;
            else
                warning('Unrecognised contour output type. Returning contDepFlags = false;');
                contDepFlags = false;
            end
            
            %If there's a problem with contour-dependence, it may be necessary
            %to flip the crack-normal vector used for contour integral
            %evaluation on one or more of the cracks
            if any(contDepFlags)
                %Show a warning and then attempt to modify the .inp file
                contDepWarnFlagLevel1(k1) = true;
                warning('Contour-dependent result detected. Attepting to flip crack-normal vector(s) in the .inp file and re-run.')
                if isfield(paramLogStruct(k1).naturalParams.outputRequests,'contourNameA') && isfield(paramLogStruct(k1).naturalParams.outputRequests,'contourNameB')
                    flipNormalStatus = int_defects_flip_normal_in_inp( [jobName,'.inp'], contDepFlags,...
                        paramLogStruct(k1).naturalParams.outputRequests.contourNameA,...
                        paramLogStruct(k1).naturalParams.outputRequests.contourNameB);
                elseif isfield(paramLogStruct(k1).naturalParams.outputRequests,'contourNameA')
                    flipNormalStatus = int_defects_flip_normal_in_inp( [jobName,'.inp'], contDepFlags,...
                        paramLogStruct(k1).naturalParams.outputRequests.contourNameA);
                elseif isfield(paramLogStruct(k1).naturalParams.outputRequests,'contourNameB')
                    flipNormalStatus = int_defects_flip_normal_in_inp( [jobName,'.inp'], contDepFlags,...
                        paramLogStruct(k1).naturalParams.outputRequests.contourNameB);
                end
                
                if flipNormalStatus
                    %Clean up Abaqus files from the first attempt (leaving the
                    %.inp file of course)
                    abaqus_cleanup;
                    
                    %Re-run the job
                    [jobExitStatus] = abaqus_run(jobName,timeOutTime,noCPUs,false,abqVer);
                    jobExitStatusSliced(k1) = jobExitStatus; %parfor sliced output variable used for returning the exitStatus of each job
                    
                    %Re-check the contour-dependence of the result
                    if strcmpi(paramLogStruct(k1).naturalParams.outputRequests.contourTypeA,'K_FACTORS')    %It is assumed that for multiple cracks, the type of contour integral output requested is the same
                        [ kCheckArray ] = read_dat_contour_integral( [jobName,'.dat'], 'k', 1);
                        [ contDepFlags ] = int_defects_contour_independence_check( kCheckArray, noContourRequests, 'k', 0.01 );
                    elseif strcmpi(paramLogStruct(k1).naturalParams.outputRequests.contourTypeA,'J_INTEGRAL')
                        [ jCheckArray ] = read_dat_contour_integral( [jobName,'.dat'], 'j', 1);
                        [ contDepFlags ] = int_defects_contour_independence_check( jCheckArray, noContourRequests, 'j', 0.02 );    %Note greater tolerance for J then for K
                    else
                        warning('Contour-dependence checking not yet implemented for T-stress results.');
                        contDepFlags = false;
                    end
                    
                    if any(contDepFlags)
                        contDepWarnFlagLevel2(k1) = true;
                        warning('Flipping the crack-normal vector did not improve contour-independence of the result... Moving on to the next model.')
                    else
                        contDepWarnFlagLevel2(k1) = false;
                        disp('Flipping the crack-normal vector appears to have improved the contour-independence of the result. Treat result for this model with caution.')
                    end
                else
                    contDepWarnFlagLevel2(k1) = true;
                    warning('Unable to flip the normal vector in the .dat file.')
                end
            else
                disp('Contour integral results for this model appear to be sufficiently path-independent.')
                contDepWarnFlagLevel1(k1) = false;
                contDepWarnFlagLevel2(k1) = false;
            end
        else
            warning('No .dat file appears to have been written for this model... Skipping.')
            contDepWarnFlagLevel1(k1) = false;
            contDepWarnFlagLevel2(k1) = false;
        end
        
        %Now that we've checked the result, you can delete all abaqus files
        %except the .dat and .sta files, if requested to.
        if abaqusCleanupFlag
            abaqus_cleanup({'.com','.inp','.log','.msg','.odb','.prt','.sim','.stt','.filepart'});
        end
        
        %Return to main directory
        cd ..
    else
        contDepWarnFlagLevel1(k1) = NaN;
        contDepWarnFlagLevel2(k1) = NaN;
    end
end
if ~contDepCheckFlag
    disp('Contour-independence checks have NOT been performed on the results.');
end

%Total time taken to run all models
modelRunTotalWalltime = toc(totalTime);
clear modelTime
clear totalTime

%% Save paramlogstruct
if length(jobExitStatusSliced) ~= length(paramLogStruct)
    error('The number of entries in jobExitStatusSliced is not the same as number of entries in paramLogStruct.')
elseif length(contDepWarnFlagLevel1) ~= length(paramLogStruct)
    error('The number of entries in contDepWarnFlagLevel1 is not the same as number of entries in paramLogStruct.')
elseif length(contDepWarnFlagLevel2) ~= length(paramLogStruct)
    error('The number of entries in contDepWarnFlagLevel2 is not the same as number of entries in paramLogStruct.')
else
    for k1 = 1:length(jobExitStatusSliced)
        paramLogStruct(k1,1).jobExitStatus = jobExitStatusSliced(k1);
        paramLogStruct(k1,1).contDepWarnFlagLevel1 = contDepWarnFlagLevel1(k1);
        paramLogStruct(k1,1).contDepWarnFlagLevel2 = contDepWarnFlagLevel2(k1);
    end
    save('workspace_dump.mat','paramLogStruct','modelRunTotalWalltime','-append') ;    %Write the ammended paramLogStruct to existing workspace_dump.mat
    disp('workspace_dump.mat has been saved.')
end

%cd back to original directory
cd ..
disp('End of int_defects_run_parametric_parallel.m');
end