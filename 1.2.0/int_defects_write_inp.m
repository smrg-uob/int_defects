function [exitStatus, varargout] = int_defects_write_inp( modelParamStruct, filenamesStruct, varargin)
%int_defects_write_inp.m
%Harry Coules 2015
%
%DESCRIPTION
%This function is used to generate an Abaqus input file using a 'master'
%Python script and a set of parameters. The parameters specified in
%naturalParamStruct are used with the master python script specified in
%filenamesStruct.filenameMasterPy to write a new python script. This is
%then executed using Abaqus/CAE kernel to produce an input (.inp) file for
%the model.
%
%INPUT ARGUMENTS
%   modelParamStruct - A structure containing crack dimensions.
%   filenamesStruct - A structure containing filenames:
%       .filenameMasterPy - name of the master python script.
%       .filenameGeneratedPy - name of the python script to be generated.
%       .filenameGeneratedInp - name of the .inp file which will be
%           generated (this will be defined by the job name in the master
%           python script).
%       .filenameInp - if specified, the .inp file will be renamed to this.
%   -OPTIONAL-
%   writeInpTimeout - The length of time (in seconds) after which the .inp
%       file write operation will be deemed to have failed.
%   abqVer - String identifying the abaqus version to use. This should
%       normally be set to 'abaqus', in which case the default installed
%       version (normally the most recent available version) will be used.
%       However, specific Abaqus versions eg. 'abq6121' or 'abq6141' can be
%       requested.
%
%OUTPUT ARGUMENTS
%   exitStatus - Logical indicating whether or not a .inp file was written
%       successfully.
%   -OPTIONAL-
%   exitMessage - A string indicating the exit status. Blank if the
%       function believed to have completed successfully. Otherwise, it
%       contains information about the reason for .inp file write failure.
%   inpWriteWarnings - A cell array of strings containing any warnings
%       returned by Abaqus CAE.
%
%NOTES
% - abqVer can only be specified if writeInpTimeout has also been
%   specified. The two optional input arguments bust be given in the order:
%   writeInpTimeout, abqVer.
%
%% Process input
%Optional input argumens
if isempty(varargin)
    writeInpTimeout = 60;   %Default to 60s
    abqVer = 'abaqus';      %Default to 'abaqus'
elseif length(varargin) == 1
    writeInpTimeout = varargin{1};
    abqVer = 'abaqus';
elseif length(varargin) == 2
    writeInpTimeout = varargin{1};
    abqVer = varargin{2};
elseif length(varargin) > 2
    error('Too many input arguments');
end

%Default Abaqus replay filename
if ~isfield(filenamesStruct,'filenameRpy');
    filenamesStruct.filenameRpy = 'abaqus.rpy';
end

%% Check for modelParamStruct.validFlag
if ~isfield(modelParamStruct,'validFlag')
    warning('Field validFlag not found in modelParamStruct. Assuming  modelParamStruct.validFlag = true.')
    modelParamStruct.validFlag = true;
end

%% Write the Python script and the .inp file
exitMessage = [];
inpWriteWarnings = {};
if modelParamStruct.validFlag
    %Write python script
    int_defects_write_py( modelParamStruct, filenamesStruct.filenameMasterPy, filenamesStruct.filenameGeneratedPy);
    
    %Check if we're running on a PC or a unix system
    if ispc
        PCflag = true;
    elseif isunix
        PCflag = false;
    else
        error('Unrecognised system type - int_defects_write_inp only works on PCs and unix systems.')
    end
    
    %Execute python script. This should write the .inp file.  
    dosCommand = [abqVer,' cae noGUI=',filenamesStruct.filenameGeneratedPy];
    if PCflag
        [~,~] = dos(dosCommand);
    elseif ~PCflag
        [~,~] = unix(dosCommand);
    else
        error('Unrecognised system type - abaqus_run only works on PCs and unix systems.')
    end
    
    %Check for timeouts
    endFlag = false;
    while ~endFlag
        %Periodically check the status file to see if the job has finished
        writeInpTime = tic;
        while toc(writeInpTime) < writeInpTimeout && ~endFlag
            pause(3);
            fid = fopen('abaqus.rpy','r','n','UTF-8');   %Try and read the .rpy file
            if fid ~= -1    %i.e. If a .rpy file has been opened
                C = textscan(fid, '%s', Inf, 'Delimiter','\n'); C = C{1};   %Read it into a cell array
                fclose(fid);
                
                %Read the lines of the .rpy file in turn
                k = 1;
                rpyScanEndFlag = false;
                inpWriteWarnings = {};  %Note: Each time, before we scan the .rpy file, clear inpWriteWarnings.
                while ~rpyScanEndFlag  %A while loop here (rather than k = 1:length(C)) means that we can stop as soon as we find one of the exit strings.
                    %Check for warnings
                    if strncmpi(char(C(k)),'#: Warning:',11)
                        inpWriteWarnings = [inpWriteWarnings;char(C(k))];
                    end
                    %Check for strings indicating write is either finished
                    %or errored-out
                    if strncmpi(char(C(k)),'#* The input file was not generated',35);  %Failed to write the .inp file
                        rpyScanEndFlag = true;
                        endFlag = true;
                        exitStatus = false;
                        exitMessage = char(C(k));
                        warning('CAE failed to write a .inp file using the .py script for this model.')
                    elseif strncmpi(char(C(k)),'#* Shell sweep feature failed',29);  %Failed on shell sweep feature - a common problem for my interacting defect models
                        rpyScanEndFlag = true;
                        endFlag = true;
                        exitStatus = false;
                        exitMessage = char(C(k));
                        warning('CAE failed to write a .inp file using the .py script for this model. nB. Failed while generating a shell sweep feature.')
                    elseif strncmpi(char(C(k)),'#* Feature creation failed',26) %General failed to create feature error
                        rpyScanEndFlag = true;
                        endFlag = true;
                        exitStatus = false;
                        exitMessage = char(C(k));
                        warning('CAE failed to write a .inp file using the .py script for this model. nB. General failure to create a feature.')
                    elseif strncmpi(char(C(k)),'#* Error:',9) %General error
                        rpyScanEndFlag = true;
                        endFlag = true;
                        exitStatus = false;
                        exitMessage = char(C(k));
                        warning('CAE failed to write a .inp file using the .py script for this model. nB. General error.')
                    elseif strncmpi(char(C(k)),'#: Warning: The following part instances are not meshed',55)
                        rpyScanEndFlag = true;
                        endFlag = true;
                        exitStatus = false;
                        exitMessage = char(C(k));
                        warning('CAE failed to write a .inp file using the .py script for this model. nB. Could not mesh.')
                    elseif strncmpi(char(C(k)),'#: Warning: The following contour integral objects contains',59)
                        rpyScanEndFlag = true;
                        endFlag = true;
                        exitStatus = false;
                        exitMessage = char(C(k));
                        warning('CAE failed to write a .inp file using the .py script for this model. nB. Problem with contour integral object.')
                    elseif strncmpi(char(C(k)),'#: RT script done',17); %Successful writing of the .inp file - note I am checking for this string last.
                        rpyScanEndFlag = true;
                        endFlag = true;
                        exitStatus = true;
                    end
                    if k >= length(C)
                        rpyScanEndFlag = true;  %Stop reading the file if you've got to the end
                    else
                        k = k+1;    %Otherwise increment the line number.
                    end
                end
            end
        end
        %Then if the write operation hasn't ended, check the timeout
        if ~endFlag
            if toc(writeInpTime) > writeInpTimeout
                %Show a warning if there is a timeout
                endFlag = true;
                exitStatus = false;
                exitMessage = ['Abaqus timeout after ',num2str(toc(writeInpTime),'%4.2f'),' seconds.'];
                warning('Timeout while attempting to write .inp file for this set of parameters.')
                %Kill remaining Abaqus process
                if PCflag
                    dos('taskkill /IM ABQcaeK.exe /F');
                    dos('taskkill /IM cmd.exe /F');
                else
                    unix('pkill -9 ABQcaeK.exe');
                end
            end
        end
    end
    
    %If required, rename the input file that has been written
    if isfield(filenamesStruct,'filenameInp')
        movefile(filenamesStruct.filenameGeneratedInp,filenamesStruct.filenameInp);
    end
else
    warning('This set of model parameters is not valid. No attempt at writing a .py or .inp file has been made.');
    exitStatus = false;
    exitMessage = 'modelParamStruct indicates this set of parameters not valid.';
end

%If requested, a message indicating the exit status can also be returned.
if nargout == 2
    varargout{1} = exitMessage;
elseif nargout == 3
    varargout{1} = exitMessage;
    varargout{2} = inpWriteWarnings;
elseif nargout > 3
    error('Too many output arguments.');
end

%Transform plate into a cylinder if required
if strcmpi(modelParamStruct.geometryType,'pipe')
    if exitStatus
        if strcmpi(modelParamStruct.pipeCrackA,'internal')
            int_defects_plate2cyl_auto(filenamesStruct.filenameGeneratedInp,'x', modelParamStruct.plateSizes.ri, true, false, modelParamStruct.plateSizes.b);
        elseif strcmpi(modelParamStruct.pipeCrackA,'external')
            int_defects_plate2cyl_auto(filenamesStruct.filenameGeneratedInp,'x', modelParamStruct.plateSizes.ri, true, true, modelParamStruct.plateSizes.b );
        else
            error('pipeCrackA should be either ''internal'' or ''external''.');
        end
    end
end

end