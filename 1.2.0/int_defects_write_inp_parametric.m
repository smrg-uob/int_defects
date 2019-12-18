function [] = int_defects_write_inp_parametric(varargin)
%int_defects_write_inp_parametric.m
%Harry Coules 2015
%
%DESCRIPTION
%This function writes Python .py files to perform a parametric analysis of
%interacting semi-elliptical cracks in a plate using Abaqus. It then uses
%Abaqus CAE to create corresponding Abaqus .inp files.
%
%USAGE
%The function can be used with no input arguments, in which case the user
%will be queried for all the information needed. Alternatively, exactly
%four input arguments can be supplied. These must be in the order:
%paramRangeFilename, filenameMasterPy, abqVer, writeInpTimeout.
%
%OPTIONAL INPUT ARGUMENTS
%  paramRangeFilename - String specifying the name of the .mat file
%       containing model parameter ranges in the structure paramRangeStruct.
%       The format of paramRangeStruct is explained the int_defects
%       documentation file.
%  filenameMasterPy - String specifying the name of the master .py file.
%  abqVer - String identifying the abaqus version to use. This should
%       normally be set to 'abaqus', in which case the default installed
%       version (normally the most recent available version) will be used.
%       However, specific Abaqus versions eg. 'abq6121' or 'abq6141' can be
%       requested.
%  writeInpTimeout - Timeout (in seconds) for the .inp file write
%       operation. Can be given as Inf if no time limit is required.
%
%OUTPUT ARGUMENTS
%   *none*
%
%NOTES
% - Although no output arguments are given, this function will create a
% directory (abaqus_param_files) which contains sub-directories for the
% individual models created, along with a dump of the MATLAB workspace and
% a logfile.
%
%% Optional input parameters, Load/ask user for general information
%Check how many input parameters have been provided. If there are zero
%input arguments, the user will be queried for the necessary information.
if nargin == 0
    manualInputFlag = true;
elseif nargin == 4
    manualInputFlag = false;
else
    error('Unexpected number of input arguments. Should be either 0 or 4.');
end

%Either get manual input, or assign variables from the input arguments.
if manualInputFlag
    %File containing parameter information
    paramRangeFilename = input('Name of the .mat file containing parameter ranges [default: paramRangeStruct.mat]:   ','s');
    if isempty(paramRangeFilename)
        paramRangeFilename = 'paramRangeStruct.mat';
    end
    
    %Filename structure
    filenamesStruct.filenameMasterPy = input('Name of the master .py file:   ','s');
    
    %Abaqus CAE version
    abqVer = input('Abaqus version to use. Should be ''abaqus'', ''abq6121'' or ''abq6141'' [default: ''abaqus'']','s');
    if isempty(abqVer)
        abqVer = 'abaqus';
    end
    
    %Timeout for .inp file write
    writeInpTimeout = input('Timeout (in seconds) for the .inp file write operation [default: 120]:   ');
    if isempty(writeInpTimeout)
        writeInpTimeout = 120;
    end
else
    %Note the input arguments must be in the correct order.
    paramRangeFilename = varargin{1};
    filenamesStruct.filenameMasterPy = varargin{2};
    abqVer = varargin{3};
    writeInpTimeout = varargin{4};
end
disp(['FILE CONTAINING PARAMETER RANGE STRUCTURE (paramRangeStruct): ',paramRangeFilename]);

load(paramRangeFilename);   %Load file containing parameter information
filenamesStruct.filenameGeneratedPy = 'IntCrackJob1_generated.py';  %Construct the rest of the filename structure. Note that these file and directory names are hard-coded.
filenamesStruct.filenameGeneratedInp = 'IntCrackJob1.inp';
filenamesStruct.mainFolderName = 'abaqus_param_files';
filenamesStruct.paramRangeFilename = paramRangeFilename;

%Type of system we're running on
if ispc
    PCflag = true;
elseif isunix
    PCflag = false;
else
    error('Unrecognised system type - int_defects_write_inp only works on PC and unix systems.')
end

%% Check geometry parameters
if ~isfield(paramRangeStruct,'geometryParams')
    error('geometryParams not defined in paramRangeStruct.')
else
    %BASIC GEOMETRIC PARAMETERS
    %Check for presence of a string indicating whether it's a plate or a pipe
    if isfield(paramRangeStruct.geometryParams,'geometryType')
        if ~any(strcmpi(paramRangeStruct.geometryParams.geometryType,{'plate','pipe'}))
            error('Invalid value of paramRangeStruct.geometryParams.geometryType. Should be ''plate'' or ''pipe''.');
        end
    else
        error('paramRangeStruct.geometryParams.geometryType not defined.');
    end
    
    %Pipe geometry - pipeCrackA is required
    if strcmpi(paramRangeStruct.geometryParams.geometryType,'pipe')
        if ~isfield(paramRangeStruct.geometryParams,'pipeCrackA')
            warning('Pipe surface on which Crack A exists (or its depth is measured from) is not defined. Assuming paramRangeStruct.geometryParams.pipeCrackA=''internal''.')
            paramRangeStruct.geometryParams.pipeCrackA = 'internal';
        end
    else
        paramRangeStruct.geometryParams.pipeCrackA = '';    %Set to empty string if geometryType is not 'pipe'
    end
    
    %Check for presence of logical indicating whether it's a single defect or a pair.
    if ~isfield(paramRangeStruct.geometryParams,'singleCrackFlag')
        warning('paramRangeStruct.geometryParams.singleCrackFlag not defined. Setting paramRangeStruct.geometryParams.singleCrackFlag = false.');
        paramRangeStruct.geometryParams.singleCrackFlag = false;
    end
    %Check that an appropriate master input file has been specified
    if paramRangeStruct.geometryParams.singleCrackFlag
        if isempty(strfind(lower(filenamesStruct.filenameMasterPy),'single'))
            warning('A single crack is being simulated but the master .py filename does not contain the string ''single''. Check master .py file.');
        end
    else
        if ~isempty(strfind(lower(filenamesStruct.filenameMasterPy),'single'))
            warning('A pair of cracks is being simulated but the master .py filename contains the string ''single''. Check master .py file.');
        end
    end
    
    %Specific checks for each case (single or pair).
    if paramRangeStruct.geometryParams.singleCrackFlag   %Single defect
        disp('This paramRangeStruct specifies parameters for single-crack models. Setting Crack B parameters to dummy values.');
        paramRangeStruct.geometryParams.aB1OveraA1 = NaN;
        paramRangeStruct.geometryParams.aB1OveraB2 = NaN;
        paramRangeStruct.geometryParams.dOverb = NaN;
        paramRangeStruct.geometryParams.oppositeFlag = false;
        paramRangeStruct.geometryParams.subsurfaceBFlag = false;
        if ~isfield(paramRangeStruct.geometryParams,'subsurfaceAFlag')
            warning('paramRangeStruct.geometryParams.subsurfaceAFlag not defined. Setting paramRangeStruct.geometryParams.subsurfaceAFlag = false.');
            paramRangeStruct.geometryParams.subsurfaceAFlag = false;
        end
        if paramRangeStruct.geometryParams.subsurfaceAFlag
            if ~isfield(paramRangeStruct.geometryParams,'d2Overa') && ~isfield(paramRangeStruct.geometryParams,'d2Overb') && ~isfield(paramRangeStruct.geometryParams,'normalisedOffsetY')
                error('Crack A is a subsurface defect, but its depth (paramRangeStruct.geometryParams.d2Overa, paramRangeStruct.geometryParams.d2Overb or paramRangeStruct.geometryParams.normalisedOffsetY) is not defined.');
            elseif sum([isfield(paramRangeStruct.geometryParams,'d2Overa'),isfield(paramRangeStruct.geometryParams,'d2Overb'),isfield(paramRangeStruct.geometryParams,'normalisedOffsetY')])>1
                error('Crack A is a subsurface defect, but its depth is defined by two or more conflicting variables (paramRangeStruct.geometryParams.d2Overa and/or paramRangeStruct.geometryParams.d2Overb and/or paramRangeStruct.geometryParams.normalisedOffsetY).');
            end
        else
            paramRangeStruct.geometryParams.d2Overb = NaN;
        end
    else   %Pairs of defects
        paramRangeStruct.geometryParams.subsurfaceAFlag = false;   %If there are multiple cracks, Crack A must be at the surface.
        if isfield(paramRangeStruct.geometryParams,'oppositeFlag')
            if isfield(paramRangeStruct.geometryParams,'subsurfaceBFlag')
                if paramRangeStruct.geometryParams.subsurfaceBFlag
                    %Check for conflicting statements regarding crack position.
                    if ~paramRangeStruct.geometryParams.oppositeFlag
                        error('Unacceptable input in paramRangeStruct: geometryParams.oppositeFlag = false but geometryParams.subsurfaceBFlag = true.')
                    end
                    %Check for presence of parameters which define the depth of
                    %embedded Crack B.
                    if ~any(isfield(paramRangeStruct.geometryParams,{'d2Overa','d2Overb','SOvera','SOverb','SOverc','SOverSqrtac','normalisedOffsetY'}))
                        error('Relative depth of subsurface crack (d2Overa, d2Overb, SOvera, SOverb, SOverc, SOverSqrtac, normalisedOffsetY) is not defined.')
                    end
                end
            else
                error('paramRangeStruct.geometryParams.subsurfaceBFlag not defined.')
            end
        else
            warning('paramRangeStruct.geometryParams.oppositeFlag is not defined. Setting to false.');
            paramRangeStruct.geometryParams.oppositeFlag = false;
            if ~isfield(paramRangeStruct.geometryParams,'subsurfaceBFlag')
                warning('paramRangeStruct.geometryParams.subsurfaceBFlag is not defined. Setting to false.');
                paramRangeStruct.geometryParams.subsurfaceBFlag = false;
            end
        end
    end
    
    %PARAMETERS DEFINING SPECIFIC ASPECTS OF THE GEOMETRY
    %Pipe radius
    if strcmpi(paramRangeStruct.geometryParams.geometryType,'pipe')
        if isfield(paramRangeStruct.geometryParams,'ri')
            paramRangeStruct.geometryParams.plDepth = paramRangeStruct.geometryParams.ri*pi;
        else
            error('Pipe internal radius (paramRangeStruct.geometryParams.ri) is not defined.');
        end
    else
        paramRangeStruct.geometryParams.ri = NaN;
    end
    
    %Wall thickness parameter (b or riOverro).
    if strcmpi(paramRangeStruct.geometryParams.geometryType,'pipe')
        if isfield(paramRangeStruct.geometryParams,'riOverro')
            if isfield(paramRangeStruct.geometryParams,'b')
                error('Only one of the parameters riOverro and b should be defined by the user.');
            else
                paramRangeStruct.geometryParams.b = (paramRangeStruct.geometryParams.ri./paramRangeStruct.geometryParams.riOverro)-paramRangeStruct.geometryParams.ri;
                k0max = length(paramRangeStruct.geometryParams.riOverro);
            end
        else
            if isfield(paramRangeStruct.geometryParams,'b')
                paramRangeStruct.geometryParams.riOverro = paramRangeStruct.geometryParams.ri./(paramRangeStruct.geometryParams.ri+paramRangeStruct.geometryParams.b);
                k0max = length(paramRangeStruct.geometryParams.b);
            else
                error('Neither riOverro nor b is defined in paramRangeStruct.geometryParams.');
            end
        end
    else
        if isfield(paramRangeStruct.geometryParams,'riOverro')
            error('For a plate model, paramRangeStruct.geometryParams.riOverro should not be defined.')
        end
        if isfield(paramRangeStruct.geometryParams,'b')
            k0max = length(paramRangeStruct.geometryParams.b);
        else
            error('For a plate model, paramRangeStruct.geometryParams.b must be defined.')
        end
    end
    
    %Variable defining the depth of Crack A (aA1Overb)
    if isfield(paramRangeStruct.geometryParams,'aA1Overb')
        k1max = length(paramRangeStruct.geometryParams.aA1Overb);
    else
        error('Crack A depth parameter paramRangeStruct.geometryParams.aA1Overb not defined.');
    end
    
    %Variable defining the aspect ratio of Crack A (aA1OveraA2)
    if isfield(paramRangeStruct.geometryParams,'aA1OveraA2')
        k2max = length(paramRangeStruct.geometryParams.aA1OveraA2);
        paramRangeStruct.geometryParams.aA1OveraA2(paramRangeStruct.geometryParams.aA1OveraA2==1)=0.9999;   %Replace Crack A aspect ration of 1 with 0.9999. Eliminates Abaqus node-numbering problems.
    else
        error('Crack A aspect ratio parameter paramRangeStruct.geometryParams.aA1OveraA2 not defined.');
    end
    
    %Check what variable is used to define the depth of Crack B in
    %paramRangeStruct.geometryParams. Get the number of parameters in this range (k2max)
    k3max = [];
    if isfield(paramRangeStruct.geometryParams,'aB1OveraA1')
        if ~isempty(k3max)
            error('paramRangeStruct.geometryParams contains multiple fields which specify the depth of Crack B.')
        else
            k3max = length(paramRangeStruct.geometryParams.aB1OveraA1);
        end
    elseif isfield(paramRangeStruct.geometryParams,'aB1Overb')
        if ~isempty(k3max)
            error('paramRangeStruct.geometryParams contains multiple fields which specify the depth of Crack B.')
        else
            k3max = length(paramRangeStruct.geometryParams.aB1Overb);
        end
    end
    
    %Variable defining the aspect ratio of Crack B (aB1OveraB2)
    if isfield(paramRangeStruct.geometryParams,'aB1OveraB2')
        k4max = length(paramRangeStruct.geometryParams.aB1OveraB2);
        paramRangeStruct.geometryParams.aB1OveraB2(paramRangeStruct.geometryParams.aB1OveraB2==1)=0.9999;   %Replace Crack B aspect ration of 1 with 0.9999. Eliminates Abaqus node-numbering problems.
    else
        error('Crack B aspect ratio parameter paramRangeStruct.geometryParams.aB1OveraB2 not defined.');
    end
    
    %Check what variable is used to define the distance between ellipses (in
    %the x-direction) in paramRangeStruct.geometryParams. Get the number of
    %parameters in this range (k5max).
    k5max = [];
    if isfield(paramRangeStruct.geometryParams,'dOverb')
        if ~isempty(k5max)
            error('paramRangeStruct.geometryParams contains multiple fields which specify the x-direction offset between cracks.')
        else
            k5max = length(paramRangeStruct.geometryParams.dOverb);
        end
    elseif isfield(paramRangeStruct.geometryParams,'dOveraA1')
        if ~isempty(k5max)
            error('paramRangeStruct.geometryParams contains multiple fields which specify the x-direction offset between cracks.')
        else
            k5max = length(paramRangeStruct.geometryParams.dOveraA1);
        end
    elseif isfield(paramRangeStruct.geometryParams,'dOverMax_aA1_2aB1')
        if ~isempty(k5max)
            error('paramRangeStruct.geometryParams contains multiple fields which specify the x-direction offset between cracks.')
        else
            k5max = length(paramRangeStruct.geometryParams.dOverMax_aA1_2aB1);
        end
    elseif isfield(paramRangeStruct.geometryParams,'twoaA2OverdCents')
        if ~isempty(k5max)
            error('paramRangeStruct.geometryParams contains multiple fields which specify the x-direction offset between cracks.')
        else
            k5max = length(paramRangeStruct.geometryParams.twoaA2OverdCents);
        end
    elseif isfield(paramRangeStruct.geometryParams,'overlapRatio')
        if ~isempty(k5max)
            error('paramRangeStruct.geometryParams contains multiple fields which specify the x-direction offset between cracks.')
        else
            k5max = length(paramRangeStruct.geometryParams.overlapRatio);
        end
    elseif isfield(paramRangeStruct.geometryParams,'normalisedOffset')
        if ~isempty(k5max)
            error('paramRangeStruct.geometryParams contains multiple fields which specify the x-direction offset between cracks.')
        else
            k5max = length(paramRangeStruct.geometryParams.normalisedOffset);
        end
    elseif isfield(paramRangeStruct.geometryParams,'sOverb')
        if ~isempty(k5max)
            error('paramRangeStruct.geometryParams contains multiple fields which specify the x-direction offset between cracks.')
        else
            k5max = length(paramRangeStruct.geometryParams.sOverb);
        end
    else
        error('paramRangeStruct.geometryParams contains neither dOverb nor dOveraA1 nor twoaA2OverdCents nor overlapRatio nor sOverb, and is not specified to be for a single crack model.')
    end
    
    %Check what variable is used to define the depth of Crack 2 in
    %paramRangeStruct. Get the number of parameters in this range (k7max).
    k7max = [];
    if isfield(paramRangeStruct.geometryParams,'d2Overa')
        if ~isempty(k7max)
            error('paramRangeStruct.geometryParams contains multiple fields which specify the depth of the subsurface crack.')
        else
            k7max = length(paramRangeStruct.geometryParams.d2Overa);
        end
    end
    if isfield(paramRangeStruct.geometryParams,'d2Overb')
        if ~isempty(k7max)
            error('paramRangeStruct.geometryParams contains multiple fields which specify the depth of the subsurface crack.')
        else
            k7max = length(paramRangeStruct.geometryParams.d2Overb);
        end
    end
    if isfield(paramRangeStruct.geometryParams,'SOvera')
        if ~isempty(k7max)
            error('paramRangeStruct.geometryParams contains multiple fields which specify the depth of the subsurface crack.')
        else
            k7max = length(paramRangeStruct.geometryParams.SOvera);
        end
    end
    if isfield(paramRangeStruct.geometryParams,'SOverb')
        if ~isempty(k7max)
            error('paramRangeStruct.geometryParams contains multiple fields which specify the depth of the subsurface crack.')
        else
            k7max = length(paramRangeStruct.geometryParams.SOverb);
        end
    end
    if isfield(paramRangeStruct.geometryParams,'SOverc')
        if ~isempty(k7max)
            error('paramRangeStruct.geometryParams contains multiple fields which specify the depth of the subsurface crack.')
        else
            k7max = length(paramRangeStruct.geometryParams.SOverc);
        end
    end
    if isfield(paramRangeStruct.geometryParams,'SOverSqrtac')
        if ~isempty(k7max)
            error('paramRangeStruct.geometryParams contains multiple fields which specify the depth of the subsurface crack.')
        else
            k7max = length(paramRangeStruct.geometryParams.SOverSqrtac);
        end
    end
    if isfield(paramRangeStruct.geometryParams,'normalisedOffsetY')
        if ~isempty(k7max)
            error('paramRangeStruct.geometryParams contains multiple fields which specify the depth of the subsurface crack.')
        else
            k7max = length(paramRangeStruct.geometryParams.normalisedOffsetY);
        end
    end
    if isempty(k7max)
        %Silently default to paramRangeStruct.geometryParams.d2Overb = NaN
        paramRangeStruct.geometryParams.d2Overb = NaN;
        k7max = length(paramRangeStruct.geometryParams.d2Overb);
    end
end

%% Check loading parameters
if ~isfield(paramRangeStruct,'loadingParams')
    error('loadingParams not defined in paramRangeStruct.')
else
    if ~isfield(paramRangeStruct.loadingParams,'loadingType')
        error('loadingParams.loadingType not defined in paramRangeStruct.');
    else
        %Valid strings for loadingType: 'load', 'moment', 'combined', 'crackFacePressure', 'pipePressure', 'pipePressureInclCrack'.
        if strcmpi(paramRangeStruct.loadingParams.loadingType,'crackFacePressure')
            if ~isfield(paramRangeStruct.loadingParams,'CFPArray')
                if isfield(paramRangeStruct.loadingParams,'CFPstring')
                    paramRangeStruct.loadingParams.CFPArray = {paramRangeStruct.loadingParams.CFPstring};
                    clear paramRangeStruct.loadingParams.CFPstring;
                else
                    error('crackFacePressure loading array/string not defined in paramRangeStruct.')
                end
            end
            k9max = size(paramRangeStruct.loadingParams.CFPArray,1);
        elseif any(strcmpi(paramRangeStruct.loadingParams.loadingType,{'load','moment','combined'}))
            if ~isfield(paramRangeStruct.loadingParams,'loadMagnitude')
                error('Magnitude of remote load not defined.');
            end
            k9max = size(paramRangeStruct.loadingParams.loadMagnitude,1);
        elseif any(strcmpi(paramRangeStruct.loadingParams.loadingType,{'pipePressure','pipePressureInclCrack'}))
            if ~isfield(paramRangeStruct.loadingParams,'pressureMagnitude')
                error('Magnitude of pipe pressure not defined.');
            end
            k9max = size(paramRangeStruct.loadingParams.pressureMagnitude,1);
        else
            error('Invalid loading mode defined in paramRangeStruct.loadingParams.loadingType.')
        end
        
        %Check biaxialFlag and biaxialStress
        if isfield(paramRangeStruct.loadingParams,'biaxialFlag')
            if paramRangeStruct.loadingParams.biaxialFlag
                if ~isfield(paramRangeStruct.loadingParams,'biaxialStress')
                    error('Magnitude of stress in ''biaxial'' (i.e. X) direction not defined.');
                end
            else
                paramRangeStruct.loadingParams.biaxialStress = 0;
            end
        else
            paramRangeStruct.loadingParams.biaxialFlag = false;
            paramRangeStruct.loadingParams.biaxialStress = 0;
        end
    end
end

%% Check material parameters
%Can define either linear elasticity, deformation plasticity or incremental plasticity.
if ~isfield(paramRangeStruct,'materialParams')
    error('materialParams not defined in paramRangeStruct.');
else
    if ~isfield(paramRangeStruct.materialParams,'type')
        error('Material type not defined in paramRangeStruct.materialParams.');
    else
        if any(strcmpi(paramRangeStruct.materialParams.type,{'linear elasticity','deformation plasticity','incremental plasticity'}))
            if isfield(paramRangeStruct.materialParams,'nu')
                k8max = length(paramRangeStruct.materialParams.nu);
            else
                error('Material Poisson''s ratio (paramRangeStruct.materialParams.nu) is undefined.');
            end
        else
            error('Material type (paramRangeStruct.materialParams.type) is unrecognised.');
        end
    end
end

%% Check modelling parameters
if ~isfield(paramRangeStruct,'modelParams')
    error('modelParams not defined in paramRangeStruct.')
else
    %Check presence of a logical indicating whether geometric nonlinearity should be considered
    if ~isfield(paramRangeStruct.modelParams,'nlgeomFlag')
        warning('paramRangeStruct.modelParams.nlgeomFlag not defined. Setting paramRangeStruct.modelParams.nlgeomFlag = false.');
        paramRangeStruct.modelParams.nlgeomFlag = false;
    end
    
    %Fields for statically-defined crack tip element sizes.
    %Note that paramRangeStruct.modelParams.rpA2 and paramRangeStruct.modelParams.rpB2 are
    %only used in int_defects_calc_model_params if paramRangeStruct.modelParams.varShellSizeFlag is both present and false.
    if ~isfield(paramRangeStruct.modelParams,'rpA2')
        paramRangeStruct.modelParams.rpA2 = NaN;
    end
    if ~isfield(paramRangeStruct.modelParams,'rpB2')
        paramRangeStruct.modelParams.rpB2 = NaN;
    end
    if length(paramRangeStruct.modelParams.rpA2) ~= length(paramRangeStruct.modelParams.rpB2)
        error('Fields rpA2 and rpB2 in paramRangeStruct.modelParams should be vectors with the same length.');
    end
    k6max = length(paramRangeStruct.modelParams.rpA2);
    
    %Check step definition string
    if ~isfield(paramRangeStruct.modelParams,'stepStr')
        paramRangeStruct.modelParams.stepStr = '';
    end
end

%% Check output requests
if ~isfield(paramRangeStruct,'outputRequests')
    error('outputRequests not defined in paramRangeStruct.');
end

%% Write parametric run files
%Make main directory and move to it
if ~logical(exist(filenamesStruct.mainFolderName,'dir'))    %Only create if it doesn't already exist
    mkdir(filenamesStruct.mainFolderName)
    if ispc
        copyfile(filenamesStruct.filenameMasterPy,['./',filenamesStruct.mainFolderName]);    %Copy the master .py file to new directory
    elseif isunix
        system(['cp ',filenamesStruct.filenameMasterPy,' ./',filenamesStruct.mainFolderName]);
    else
        error('Unrecognised system type')
    end
else   %Otherwise show a warning
    warning(['Folder ',filenamesStruct.mainFolderName,' already exists. Using existing folder rather than creating a new one.']);
end
cd(filenamesStruct.mainFolderName);

%Step through the parameter space and write .inp files to individual directories
k = 0;
kTotal = k0max*k1max*k2max*k3max*k4max*k5max*k6max*k7max*k8max*k9max;
totalTime = tic;
disp(' ');
for k0 = 1:k0max
    for k1 = 1:k1max
        for k2 = 1:k2max
            for k3 = 1:k3max
                for k4 = 1:k4max
                    for k5 = 1:k5max
                        for k6 = 1:k6max
                            for k7 = 1:k7max
                                for k8 = 1:k8max
                                    for k9 = 1:k9max
                                        %Run number k.
                                        modelTime = tic;
                                        k = k+1;
                                        
                                        %Create directory and cd to it
                                        if ~logical(exist(num2str(k,'%06i'),'dir'))    %Only create if it doesn't already exist
                                            mkdir(num2str(k,'%06i'))
                                            if ispc
                                                copyfile(filenamesStruct.filenameMasterPy,['./',num2str(k,'%06i')]);    %Copy file to new directory
                                            elseif isunix
                                                system(['cp ',filenamesStruct.filenameMasterPy,' ./',num2str(k,'%06i')]);
                                            else
                                                error('Unrecognised system type')
                                            end
                                        else
                                            disp(['Folder ',num2str(k,'%06i'),' already exists.'])
                                        end
                                        cd(num2str(k,'%06i'));
                                        
                                        %Write input file into directory...
                                        if isempty(dir('*.inp'))    %... only do this if there isn't one in there already.
                                            oldInpPresent = false;
                                            disp(datestr(clock));
                                            disp(['WRITING .PY AND .INP FILES FOR PARAMETRIC MODEL ',num2str(k),' OF ',num2str(kTotal)]);
                                            
                                            %CONSTRUCT PARAMETER STRUCTURE FOR THIS RUN
                                            naturalParamStruct = paramRangeStruct;
                                            
                                            %Wall thickness or pipe aspect ratio
                                            if isfield(paramRangeStruct.geometryParams,'b')
                                                naturalParamStruct.geometryParams.b = paramRangeStruct.geometryParams.b(k0);                      %Wall thickness
                                            end
                                            if isfield(paramRangeStruct.geometryParams,'riOverro')
                                                naturalParamStruct.geometryParams.riOverro = paramRangeStruct.geometryParams.riOverro(k0);        %Pipe thickness ratio
                                            end
                                            
                                            naturalParamStruct.geometryParams.aA1Overb = naturalParamStruct.geometryParams.aA1Overb(k1);
                                            naturalParamStruct.geometryParams.aA1OveraA2 = naturalParamStruct.geometryParams.aA1OveraA2(k2);
                                            
                                            %Size of Crack B (in through-thickness dir.) - can be defined by different variables
                                            if isfield(naturalParamStruct.geometryParams,'aB1OveraA1')
                                                naturalParamStruct.geometryParams.aB1OveraA1 = naturalParamStruct.geometryParams.aB1OveraA1(k3);
                                            end
                                            if isfield(naturalParamStruct.geometryParams,'aB1Overb')
                                                naturalParamStruct.geometryParams.aB1Overb = naturalParamStruct.geometryParams.aB1Overb(k3);
                                            end
                                            
                                            %Crack B aspect ratio
                                            if isnan(naturalParamStruct.geometryParams.aB1OveraB2(k4))     %If aB1OveraB2 is set as NaN, set it the same as aA1OveraA2
                                                naturalParamStruct.geometryParams.aB1OveraB2 = naturalParamStruct.geometryParams.aA1OveraA2;
                                            else
                                                naturalParamStruct.geometryParams.aB1OveraB2 = naturalParamStruct.geometryParams.aB1OveraB2(k4);
                                            end
                                            
                                            %Distance between cracks along the plate - can be defined by different variables
                                            if isfield(naturalParamStruct.geometryParams,'dOverb')
                                                naturalParamStruct.geometryParams.dOverb = naturalParamStruct.geometryParams.dOverb(k5);
                                            end
                                            if isfield(naturalParamStruct.geometryParams,'dOveraA1')
                                                naturalParamStruct.geometryParams.dOveraA1 = naturalParamStruct.geometryParams.dOveraA1(k5);
                                            end
                                            if isfield(naturalParamStruct.geometryParams,'dOverMax_aA1_2aB1')
                                                naturalParamStruct.geometryParams.dOverMax_aA1_2aB1 = naturalParamStruct.geometryParams.dOverMax_aA1_2aB1(k5);
                                            end
                                            if isfield(naturalParamStruct.geometryParams,'twoaA2OverdCents')
                                                naturalParamStruct.geometryParams.twoaA2OverdCents = naturalParamStruct.geometryParams.twoaA2OverdCents(k5);
                                            end
                                            if isfield(naturalParamStruct.geometryParams,'overlapRatio')
                                                naturalParamStruct.geometryParams.overlapRatio = naturalParamStruct.geometryParams.overlapRatio(k5);
                                            end
                                            if isfield(naturalParamStruct.geometryParams,'sOverb')
                                                naturalParamStruct.geometryParams.sOverb = naturalParamStruct.geometryParams.sOverb(k5);
                                            end
                                            if isfield(naturalParamStruct.geometryParams,'normalisedOffset')
                                                naturalParamStruct.geometryParams.normalisedOffset = naturalParamStruct.geometryParams.normalisedOffset(k5);
                                            end
                                            
                                            %Crack tip mesh zone size
                                            naturalParamStruct.modelParams.rpA2 = naturalParamStruct.modelParams.rpA2(k6);
                                            naturalParamStruct.modelParams.rpB2 = naturalParamStruct.modelParams.rpB2(k6);
                                            
                                            %Depth of Crack B - can be defined by different variables
                                            if isfield(naturalParamStruct.geometryParams,'d2Overa')
                                                naturalParamStruct.geometryParams.d2Overa = naturalParamStruct.geometryParams.d2Overa(k7);
                                            end
                                            if isfield(naturalParamStruct.geometryParams,'d2Overb')
                                                naturalParamStruct.geometryParams.d2Overb = naturalParamStruct.geometryParams.d2Overb(k7);
                                            end
                                            if isfield(naturalParamStruct.geometryParams,'SOvera')
                                                naturalParamStruct.geometryParams.SOvera = naturalParamStruct.geometryParams.SOvera(k7);
                                            end
                                            if isfield(naturalParamStruct.geometryParams,'SOverb')
                                                naturalParamStruct.geometryParams.SOverb = naturalParamStruct.geometryParams.SOverb(k7);
                                            end
                                            if isfield(naturalParamStruct.geometryParams,'SOverc')
                                                naturalParamStruct.geometryParams.SOverc = naturalParamStruct.geometryParams.SOverc(k7);
                                            end
                                            if isfield(naturalParamStruct.geometryParams,'SOverSqrtac')
                                                naturalParamStruct.geometryParams.SOverSqrtac = naturalParamStruct.geometryParams.SOverSqrtac(k7);
                                            end
                                            if isfield(naturalParamStruct.geometryParams,'normalisedOffsetY')
                                                naturalParamStruct.geometryParams.normalisedOffsetY = naturalParamStruct.geometryParams.normalisedOffsetY(k7);
                                            end
                                            
                                            %Poisson's ratio
                                            naturalParamStruct.materialParams.nu = naturalParamStruct.materialParams.nu(k8);
                                            
                                            %Loading parameters
                                            if strcmpi(naturalParamStruct.loadingParams.loadingType,'crackFacePressure')
                                                naturalParamStruct.loadingParams.CFPstring = naturalParamStruct.loadingParams.CFPArray{k9};
                                            elseif any(strcmpi(naturalParamStruct.loadingParams.loadingType,{'load','pressure','combined'}))
                                                naturalParamStruct.loadingParams.loadMagnitude = naturalParamStruct.loadingParams.loadMagnitude(k9,:);
                                            elseif any(strcmpi(naturalParamStruct.loadingParams.loadingType,{'pipePressure','pipePressureInclCrack'}))
                                                naturalParamStruct.loadingParams.pressureMagnitude = naturalParamStruct.loadingParams.pressureMagnitude(k9,:);
                                            else
                                                error('Unrecognised loading type (paramRangeStruct.loadingParams.loadingType).');
                                            end
                                            
                                            %WRITE MODEL INPUT FILES                                            
                                            %Calculate model parameters from natural parameters
                                            [ modelParamStruct, modelParamsValidityComment ] = int_defects_calc_model_params( naturalParamStruct );
                                            
                                            %Write .py and .inp files
                                            int_defects_write_py( modelParamStruct, filenamesStruct.filenameMasterPy, filenamesStruct.filenameGeneratedPy);
                                            [inpStatus, inpWriteExitMessage, inpWriteWarnings] = int_defects_write_inp( modelParamStruct, filenamesStruct, writeInpTimeout, abqVer);
                                            
                                            %Kill remaining Abaqus process
                                            if PCflag
                                                [~,~] = dos('taskkill /IM ABQcaeK.exe /F');
                                                [~,~] = dos('taskkill /IM cmd.exe /F');
                                            else
                                                [~,~] = unix('pkill -9 ABQcaeK.exe');
                                            end
                                            
                                            %If the model parameters are valid but CAE was
                                            %not able to write an input file, there may
                                            %have been a meshing problem. Try relaxing the
                                            %mesh constraints.
                                            if modelParamStruct.validFlag && ~inpStatus
                                                %Show a warning and set modelParamStruct.relaxMeshFlag1
                                                warning('There may have been a mesh generation problem. Relaxing mesh constraints (level 1) and retrying. This may create an overly refined mesh.');
                                                modelParamStruct.relaxMeshFlag1 = true;
                                                
                                                %Delete any existing .rec and .rpy files
                                                warning('off');
                                                delete('*.rpy');
                                                delete('*.rec');
                                                warning('on');
                                                
                                                %Write .py and .inp files
                                                int_defects_write_py( modelParamStruct, filenamesStruct.filenameMasterPy, filenamesStruct.filenameGeneratedPy);
                                                [inpStatus, inpWriteExitMessage, inpWriteWarnings] = int_defects_write_inp( modelParamStruct, filenamesStruct, writeInpTimeout, abqVer);
                                                
                                                %Kill remaining Abaqus process
                                                disp('Terminating remaining Abaqus processes...');
                                                if PCflag
                                                    [~,~] = dos('taskkill /IM ABQcaeK.exe /F');
                                                    [~,~] = dos('taskkill /IM cmd.exe /F');
                                                else
                                                    [~,~] = unix('pkill -9 ABQcaeK.exe');
                                                end
                                                
                                                %Warning message to indicate input file status
                                                if inpStatus
                                                    warning('Succeeded in writing .inp after relaxing mesh constraints (level 1). Results should be treated with caution.')
                                                else
                                                    warning('Failed to write .inp after relaxing mesh constraints (level 1). Relaxing mesh constraints (level 2) and retrying. This may create an overly refined mesh.')
                                                    modelParamStruct.relaxMeshFlag2 = true;
                                                    
                                                    %Delete any existing .rec and .rpy files
                                                    warning('off')
                                                    delete('*.rpy');
                                                    delete('*.rec');
                                                    warning('on')
                                                    %Write .py and .inp files
                                                    int_defects_write_py( modelParamStruct, filenamesStruct.filenameMasterPy, filenamesStruct.filenameGeneratedPy);
                                                    [inpStatus, inpWriteExitMessage, inpWriteWarnings] = int_defects_write_inp( modelParamStruct, filenamesStruct, writeInpTimeout, abqVer);
                                                    
                                                    %Kill remaining Abaqus process
                                                    disp('Terminating remaining Abaqus processes...');
                                                    if PCflag
                                                        [~,~] = dos('taskkill /IM ABQcaeK.exe /F');
                                                        [~,~] = dos('taskkill /IM cmd.exe /F');
                                                    else
                                                        [~,~] = unix('pkill -9 ABQcaeK.exe');
                                                    end
                                                    
                                                    %Warning message to indicate input file status
                                                    if inpStatus
                                                        warning('Succeeded in writing .inp after relaxing mesh constraints (level 2). Results should be treated with caution.');
                                                    else
                                                        warning('Failed to write .inp after relaxing mesh constraints (level 2).');
                                                    end
                                                end
                                            end
                                            
                                            %Delete master file
                                            if ~isempty(dir(filenamesStruct.filenameMasterPy))
                                                delete(filenamesStruct.filenameMasterPy);
                                            end
                                            
                                            %Save parameters for this run in a structure containing
                                            %parameters for all runs
                                            paramLogStruct(k,1).filenames = filenamesStruct;
                                            paramLogStruct(k,1).modelParams = modelParamStruct;
                                            paramLogStruct(k,1).naturalParams = naturalParamStruct;
                                            paramLogStruct(k,1).inpStatus = inpStatus;
                                            paramLogStruct(k,1).modelParamsValidityComment = modelParamsValidityComment;
                                            paramLogStruct(k,1).inpWriteExitMessage = inpWriteExitMessage;
                                            paramLogStruct(k,1).inpWriteWarnings = inpWriteWarnings;
                                        else
                                            oldInpPresent = true;
                                            disp(['There is already a .inp file in folder: ',num2str(k,'%06i'),'. Not attempting to write a new .inp file for this .']);
                                        end
                                        
                                        %Modify .inp files with custom output requests if required
                                        if isfield(paramRangeStruct.outputRequests,'custom')
                                            if inpStatus
                                                find_replace_string( '** OUTPUT REQUESTS', paramRangeStruct.outputRequests.custom, 'IntCrackJob1.inp' );
                                            end
                                        end
                                        
                                        %cd back to main directory
                                        cd ..
                                        
                                        %Display
                                        if oldInpPresent
                                            disp('.INP FILE FOR THIS MODEL ALREADY EXISTED - no attempt made to modify it.');
                                        else
                                            if inpStatus
                                                disp('.INP FILE WRITE OPERATION COMPLETE - the input file was written successfully.');
                                            else
                                                disp('.INP FILE WRITE OPERATION FAILED - the input file has not been written, or may be incomplete.');
                                            end
                                        end
                                        disp(['Walltime for this model: ',num2str(toc(modelTime),'%4.2f'),' seconds']);
                                        disp(['Total walltime: ',num2str(toc(totalTime),'%4.2f'),' seconds']);
                                        disp(' ');
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%Total time taken to create all models
modelCreateTotalWalltime = toc(totalTime);
clear modelTime
clear totalTime

%Save workspace, delete master file
if ~isempty(dir('workspace_dump.mat'))
    warning('workspace_dump.mat already exists. Renaming the existing file to workspace_dump_0.mat')
    movefile('workspace_dump.mat','workspace_dump_0.mat')
end
save('workspace_dump.mat');
delete(filenamesStruct.filenameMasterPy);

%Write logfile
int_defects_write_logfile_1( paramLogStruct );

%cd back to original directory
cd ..
end