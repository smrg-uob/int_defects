function [varargout] = abaqus_run(jobName,timeOutTime,varargin)
%abaqus_run.m
%Harry Coules 2013
%
%DESCRIPTION
%This function provides a more robust method of running Abaqus analyses
%from Matlab, allowing the user to set a timeout which kills and restarts
%the job if it is taking too long.
%
%INPUT ARGUMENTS
% jobName - Name of the Abaqus job
% timeoutTime - Timeout limit in seconds. Can be set as Inf if necessary.
%
%  *OPTIONAL INPUT ARGUMENTS*
% forFileName - String defining a user subroutine fortran file, if the
%               Abaqus model is intended to be executed using a user
%               subroutine.
% noCPUs - The number of CPUs to use, if running on a
%               multi-core system. If noCPUs is excluded or given as NaN,
%               the number of CPUs will not be specified in the command
%               used to run the job.
% allowRetries - Logical to specify whether the script should try to
%               resubmit the job in case of a timeout. Default is 'true' if
%               this input argument is not supplied.
% abqVer - Abaqus version to use. Should be either 'abaqus' if the system's
%               default version (usually the most recent) should be used,
%               or a specific version eg. 'abq6121' or 'abq6141'. Default
%               is 'abaqus'.
%
%OUTPUT ARGUMENTS
%   MANDATORY:
% *none*
%
%   *OPTIONAL OUTPUT ARGUMENTS*
%exitStatus - Numeric indicating whether or not the analysis as intended.
%   0 = Did not complete
%   1 = Completed successfully
%
%NOTES
% - If more than one optional input argument is given, they must be
%   in the order: forFileName, noCPUs, allowRetries, abqVer (any may be
%   omitted from this list).
%
%% Check input arguments
%Check number of output arguments is okay before running the script
if nargout > 1
    error('Too many output arguments.')
end

allowRetries = true;
if nargin < 2
    error('Not enough input arguments')
elseif nargin == 3
    if ischar(varargin{1})
        forFileName = varargin{1};
    elseif isnumeric(varargin{1})
        noCPUs = varargin{1};
    elseif islogical(varargin{1})
        allowRetries = varargin{1};
    else
        error('When one optional input argument is provided it must be either char (forFileName), numeric (noCPUs), or logical (allowRetries).');
    end
    abqVer = 'abaqus';
elseif nargin == 4
    if ischar(varargin{1})
        forFileName = varargin{1};
        if isnumeric(varargin{2})
            noCPUs = varargin{2};
        elseif islogical(varargin{2})
            allowRetries = varargin{2};
        else
            error('When two optional input arguments are given, the second should be either numeric (for noCPUs) or logical (for allowRetries).');
        end
    elseif isnumeric(varargin{1})
        noCPUs = varargin{1};
        if islogical(varargin{2})
            allowRetries = varargin{2};
        else
            error('When two optional input arguments are given and the first is noCPUs, the second must be a logical (for allowRetries).');
        end
    else
        error('When two optional input arguments are given, the first should be either char (for forFileName) or numeric (for noCPUs).'); 
    end
    abqVer = 'abaqus';
elseif nargin == 5
    if ischar(varargin{1})
        forFileName = varargin{1};
        noCPUs = varargin{2};
        allowRetries = varargin{3};
        abqVer = 'abaqus';
    elseif isnumeric(varargin{1})
        noCPUs = varargin{1};
        allowRetries = varargin{2};
        abqVer = varargin{3};
    else
        error('When three optional input arguments are given, the first should be either char (for forFileName), or numeric (for noCPUs).');
    end
elseif nargin == 6
    if ischar(varargin{1})
        forFileName = varargin{1};
    else
        error('When four optional input arguments are given, the first must be a char (for forFileName).');
    end
    if isnumeric(varargin{2})
        noCPUs = varargin{2};
    else
        error('When four optional input arguments are given, the second must be numeric (for noCPUs).');
    end
    if islogical(varargin{3})
        allowRetries = varargin{3};
    else
        error('When four optional input arguments are given, the third must be a logical (for allowRetries).');
    end
    if ischar(varargin{4})
        abqVer = varargin{4};
    else
        error('When four optional input arguments are given, the fourth must be a char (for abqVer).');
    end
else
    error('Too many input arguments.')
end

%Check Abaqus version is valid
if ~strcmp(abqVer,'abaqus') && ~strcmp(abqVer,'abq6121') && ~strcmp(abqVer,'abq6141')
    warning('abqVer should be ''abaqus'', ''abq6121'' or ''abq6141''.');
end

%If noCPUs is given as NaN, then delete this variable so that the dos
%command string will not contain 'cpus='
if ~exist('noCPUs','var')
    noCPUs = NaN;
end

%Check if we're running on a PC or a unix system
if ispc
    PCflag = true;
elseif isunix
    PCflag = false;
else
    error('Unrecognised system type - abaqus_run only works on PCs and unix systems.')
end

%Compile the dos command to be run as a string
if exist('forFileName','var') && ~isnan(noCPUs)
    DosCommand = strcat(abqVer,' analysis job=',jobName,' user=',forFileName,' cpus=',num2str(noCPUs),' scratch="scratch" interactive');   %With user subroutine and number of CPUs
elseif ~isnan(noCPUs)
    DosCommand = strcat(abqVer,' analysis job=',jobName,' cpus=',num2str(noCPUs),' scratch="scratch" interactive');   %With number of CPUs only
elseif exist('forFileName','var')
    DosCommand = strcat(abqVer,' analysis job=',jobName,' user=',forFileName,' scratch="scratch" interactive');   %With user subroutine only
else
    DosCommand = strcat(abqVer,' analysis job=',jobName,' scratch="scratch" interactive');    %Without user subroutine of number of CPUs
end

DosCommandTerm = strcat(abqVer,' terminate job=',jobName);
filenameDat = strcat(jobName,'.dat');
filenameSta = strcat(jobName,'.sta');

%% Create scratch directory and run analysis
mkdir('scratch');
if PCflag
    [~,~] = dos(DosCommand,'-echo');
elseif ~PCflag
    [~,~] = unix(DosCommand,'-echo');
else
    error('Unrecognised system type - abaqus_run only works on PCs and unix systems.')
end

%Check for timeouts
timeOuts = 0;
endFlag = false;
exitStatus = 0; %Numeric so that in future it could be used to indicate exit statuses other than completed/did not complete.
while endFlag == false
    %Periodically check the status file to see if the job has finished
    time = 0;
    tic
    while time < timeOutTime && endFlag == false
        time = toc; %elapsed computation walltime
        pause(1);

        fid = fopen(filenameSta,'r','n','UTF-8');   %Try and read the status file
        if fid ~= -1    %i.e. If a .sta file has been opened
            C = textscan(fid, '%s', Inf, 'Delimiter','\n'); C = C{1};   %Read it into a cell array
            fclose(fid);
            for k = 1:length(C)
                if strcmpi(char(C(k)),'THE ANALYSIS HAS COMPLETED SUCCESSFULLY');
                    endFlag = true;
                    exitStatus = 1;
                elseif strcmpi(char(C(k)),'THE ANALYSIS HAS NOT BEEN COMPLETED');
                    endFlag = true;
                    exitStatus = 0;
                    warning(['Analysis of the Abaqus job ',jobName,' was not completed successfully.']);    %This allows abaqus_run to complete even if Abaqus cannot run the job
                end
            end
            
        else
            fid2 = fopen(filenameDat,'r','n','UTF-8');   %If you can't find the .sta file, look for the .dat file
            if fid2 ~= -1    %i.e. If a .dat file has been opened
                C2 = textscan(fid2, '%s', Inf, 'Delimiter','\n'); C2 = C2{1};   %Read it into a cell array
                fclose(fid2);
                for k = 1:length(C2)
                    if strcmpi(char(C2(k)),'ANALYSIS COMPLETE');
                        endFlag = true;
                        exitStatus = 1;
                    end
                end
            end
        end
    end
    
    %What to do in case of a timeout
    if time >= timeOutTime
        timeOuts = timeOuts + 1;
        if ~allowRetries
            endFlag = true;
            exitStatus = 0;
            warning('Abaqus job timeout. Retry will not be attempted because allowRetries=false. Aborting script.')
        elseif timeOuts >=6
            endFlag = true;
            exitStatus = 0;
            warning('Abaqus job timeout. Fifth retry failed. Aborting script.')
        else
            warning('Abaqus job timeout. Terminating the job, will then retry submitting.');
            
            %Try to terminate Abaqus job
            if PCflag
                dos(DosCommandTerm);
            else
                unix(DosCommandTerm);
            end
            pause(30)
            
            %Kill any remaining Abaqus processes
            if PCflag
                dos('taskkill /IM pre.exe /F');
                pause(30)
                dos('taskkill /IM abq6121.exe /F');
                pause(30)
            else
                warning('This script does not attempt to kill individual processes when running under unix.')
            end
                        
            abaqus_cleanup(jobName);    %Clean up - note that the .for file is left alone
            
            pause(30);
            
            warning('Job successfully terminated. Resubmitting the job.')
            mkdir('scratch');   %Retry
            if PCflag
                dos(DosCommand);
            else
                unix(DosCommand);
            end
        end
    end
end

%% Output exitStatus if required
if nargout == 1
    varargout{1} = exitStatus;
elseif nargout > 1
    error('Too many output arguments.');
end

end