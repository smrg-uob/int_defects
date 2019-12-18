function int_defects_calc_interaction_factor( mainFolderName, singleFolderName, varargin)
%int_defects_calc_interaction_factor.m
%Harry Coules 2015
%
%DESCRIPTION
%This function calculates the "interaction factor" for interacting
%semi-elliptical surface-breaking cracks. That is, the ratio of the Mode I
%SIF for a crack in the presence of a second crack to the Mode I SIF for
%the first crack alone. It saves the results to a .mat file.
%
%INPUT ARGUMENTS
%   mainFolderName - Name of the folder containing the data for models of
%       interacting pairs of semi-elliptical cracks.
%   singleFolderName - Name of the folder containing the data for the
%       models of single semi-elliptical cracks, as a string.
%OPTIONAL INPUT ARGUMENT
%   singleFolderName2 - A second folder containing results for single
%       elliptical subsurface cracks, to be used when investigating
%       surface-subsurface defect interactions.
%
%OUTPUT ARGUMENTS
%   *NONE* (Although a .mat file containing the results will be saved to
%   the current directory).
%
%NOTES
% - All models of interacting defects must have exactly two defects
%   in them. However, there are two methods for determining the SIF
%   distributions on solitary defects that I have used:
%       1. A similar model to the one for proximate cracks is used, but
%       the two defects are separated by a large distance. This has the
%       advantage of being simple to implement, but only if one such model
%       is run for every pair of proximate cracks - which is very
%       inefficient.
%       2. Models containing a single crack are used. When this method is
%       used, it is necessary to match up crack geometries when calculating
%       the interaction factor results in this script. In this script, the
%       logical "singleCrackFlag" is used to indicate whether the results
%       for single cracks are for the "single crack model" or "crack pair 
%       model" method, and this is determined automatically.
% - When Method 2 above is used for models containing one surface and one
%   embedded flaw, or for two surface flaws on opposite surfaces of a pipe
%   two sets of single crack models are required (i.e. one for surface and
%   one for embedded, or one for interior and one for exterior) so the
%   optional input parameter singleFolderName2 must be provided.
% - For some numeric parameters, a tolerance is used in equality checks.
%   This tolerance is hard-coded to 1 part in 1000.
%
%% Check the optional input argument
% The string singleFolderName2 is used when there are different output
%folders for comparison with Crack A and Crack B. This is used when Crack B
%is subsurface and A is not, or if Crack B is on the opposise surface of a
%pipe.
if isempty(varargin)
    twoSingleFoldersFlag = false;   %Logical to indicate whether there are two sets of single-crack results
elseif length(varargin)==1
    singleFolderName2 = varargin{1};
    twoSingleFoldersFlag = true;
else
    error('Too many input arguments.');
end

%% Load FEA results
%Load results from first set of models of single cracks
disp(' ')
disp('Loading results for remote cracks...')
cd(singleFolderName);
singleStruct = load('workspace_dump.mat', 'paramLogStruct', 'paramRangeStruct', 'outputArray');
if ~isfield(singleStruct,'paramLogStruct')
    error('Unable to load paramLogStruct from workspace_dump.mat for the set of results for a single crack.')
elseif ~isfield(singleStruct,'outputArray')
    error('Unable to load outputArray from workspace_dump.mat for the set of results for a single crack. This may be because int_defects_read_parametric (or int_defects_read_user) has not been run for this model set.')
end
if length(singleStruct.paramLogStruct) ~= length(singleStruct.outputArray)
    error('singleStruct.paramLogStruct and singleStruct.outputArray have different lengths.')
end
if isfield(singleStruct.paramLogStruct(1).naturalParams.geometryParams,'singleCrackFlag')  %Check to see if the models of remote cracks are single-crack models
    if singleStruct.paramLogStruct(1).naturalParams.geometryParams.singleCrackFlag
        singleCrackFlag = true;
    else
        singleCrackFlag = false;
    end
else
    error('Missing parameter singleCrackFlag in singleStruct.paramLogStruct(1).naturalParams.geometryParams.');
end
cd ..

%Load results from additional set of models of single cracks, if two sets of single-crack models have been used
if twoSingleFoldersFlag
    if singleCrackFlag == true
        disp('Second directory of remote crack results specified. Loading results from this as well...')
        cd(singleFolderName2);
        singleStruct2 = load('workspace_dump.mat', 'paramLogStruct', 'paramRangeStruct', 'outputArray');
        if ~isfield(singleStruct2,'paramLogStruct')
            error('Unable to load paramLogStruct from workspace_dump.mat for the set of results for a single crack.')
        elseif ~isfield(singleStruct2,'outputArray')
            error('Unable to load outputArray from workspace_dump.mat for the set of results for a single crack. This may be because int_defects_read_parametric (or int_defects_read_user) has not been run for this model set.')
        end
        if length(singleStruct2.paramLogStruct) ~= length(singleStruct2.outputArray)
            error('singleStruct2.paramLogStruct and singleStruct2.outputArray have different lengths.')
        end
    else
        error('Two folders for single crack models specified, but singleCrackFlag = false.')
    end
    cd ..
end

disp('DONE!')

%Load results for models of pairs of cracks
disp('Loading results for pairs of cracks...')
cd(mainFolderName);
mainStruct = load('workspace_dump.mat', 'paramLogStruct', 'paramRangeStruct', 'outputArray');
if ~isfield(mainStruct,'paramLogStruct')
    error('Unable to load paramLogStruct from workspace_dump.mat for main set of results.')
elseif ~isfield(mainStruct,'outputArray')
    error('Unable to load outputArray from workspace_dump.mat for main set of results. This may be because int_defects_read_parametric (or int_defects_read_user) has not been run for this model set.')
end
if length(mainStruct.paramLogStruct) ~= length(mainStruct.outputArray)
    error('mainStruct.paramLogStruct and mainStruct.outputArray have different lengths.')
end
cd ..
disp('DONE!')

%% Determine relevant single crack result and calculate interaction factors
disp('Calculating interaction factors...')
if singleCrackFlag
    for k1 = 1:length(mainStruct.outputArray)
        %Display the model number for every 100th model
        if mod(k1,100) == 0
            disp(['Current model no.: ',num2str(k1)])
        end

        %Dynamically unpack mainStruct.paramLogStruct(k1,1).naturalParams
        m = struct;
        fieldNames = fieldnames(mainStruct.paramLogStruct(k1,1).naturalParams);
        for j1 = 1:length(fieldNames)
            if isstruct(eval(['mainStruct.paramLogStruct(k1,1).naturalParams.',fieldNames{j1},';']))
                subFieldNames = fieldnames(eval(['mainStruct.paramLogStruct(k1,1).naturalParams.',fieldNames{j1},';']));
                for j2 = 1:length(subFieldNames)
                    if isfield(m,subFieldNames{j2})
                        error(['Error unpacking mainStruct.paramLogStruct(k1,1).naturalParams. Duplicate field name: ',subFieldNames{j2}]);
                    else
                        eval(['m.',subFieldNames{j2},' = mainStruct.paramLogStruct(k1,1).naturalParams.',fieldNames{j1},'.',subFieldNames{j2},';']);
                    end
                end
            else
                if isfield(m,fieldNames{j1})
                    error(['Error unpacking mainStruct.paramLogStruct(k1,1).naturalParams. Duplicate field name: ',fieldNames{j1}]);
                else
                    eval(['m.',fieldNames{j1},' = mainStruct.paramLogStruct(k1,1).naturalParams.',fieldNames{j1},';']);
                end
            end
        end
        
        %Work out which models for a single crack to use from the naturalParams
        mainStruct.paramLogStruct(k1,1).singleCrackNo = [0,0];  %nB. This is in the format: [no. of single crack corresponding to Crack A, no. of single crack corresponding to Crack B]
        
        %For Crack A
        k2 = 1;
        while mainStruct.paramLogStruct(k1,1).singleCrackNo(1) == 0
            %Dynamically unpack singleStruct.paramLogStruct(k2,1).naturalParams
            s = struct;
            fieldNames = fieldnames(singleStruct.paramLogStruct(k2,1).naturalParams);
            for j1 = 1:length(fieldNames)
                if isstruct(eval(['singleStruct.paramLogStruct(k2,1).naturalParams.',fieldNames{j1},';']))
                    subFieldNames = fieldnames(eval(['singleStruct.paramLogStruct(k2,1).naturalParams.',fieldNames{j1},';']));
                    for j2 = 1:length(subFieldNames)
                        if isfield(s,subFieldNames{j2})
                            error(['Error unpacking singleStruct.paramLogStruct(k2,1).naturalParams. Duplicate field name: ',subFieldNames{j2}]);
                        else
                            eval(['s.',subFieldNames{j2},' = singleStruct.paramLogStruct(k2,1).naturalParams.',fieldNames{j1},'.',subFieldNames{j2},';']);
                        end
                    end
                else
                    if isfield(s,fieldNames{j1})
                        error(['Error unpacking singleStruct.paramLogStruct(k2,1).naturalParams. Duplicate field name: ',fieldNames{j1}]);
                    else
                        eval(['s.',fieldNames{j1},' = singleStruct.paramLogStruct(k2,1).naturalParams.',fieldNames{j1},';']);
                    end
                end
            end
            
            %Check for model equivalence
            equalFlagGeom = strcmpi(m.geometryType,s.geometryType);             %General Geometry Type
            if strcmpi(m.geometryType,'pipe')
                equalFlagGeom = equalFlagGeom && strcmpi(m.pipeCrackA,s.pipeCrackA);	%Crack A should be on same pipe surface.
                equalFlag0 = isequal(m.ri,s.ri) && isequal(m.b,s.b);                    %Pipe radius & wall thickness
            else
                equalFlag0 = isequal(m.b,s.b);                                  %Plate thickness
            end
            equalFlag1 = isequal(m.aA1Overb,s.aA1Overb);                        %Crack depth
            equalFlag2 = abs(m.aA1OveraA2-s.aA1OveraA2) < 1e-3*s.aA1OveraA2;    %Crack aspect ratio
            equalFlag3 = isequal(m.nu,s.nu);                                    %Poisson's ratio
            if strcmpi(m.loadingType,'crackFacePressure')                       %Loading conditions
                equalFlag4 = strcmpi(m.CFPstring,s.CFPstring);
            elseif any(strcmpi(m.loadingType,{'pipePressure','pipePressureInclCrack'}))
                equalFlag4 = strcmpi(m.loadingType,s.loadingType) && isequal(m.pressureMagnitude,s.pressureMagnitude);
            else
                equalFlag4 = strcmpi(m.loadingType,s.loadingType) && isequal(m.loadMagnitude,s.loadMagnitude);
            end
            equalFlag5 = isequal(m.biaxialFlag,s.biaxialFlag) && (isequal(m.biaxialStress,s.biaxialStress) || (isnan(m.biaxialStress) && isnan(s.biaxialStress)));  %Loading conditions - biaxial direction
            %Check that all parameters are (simultaneously) correct
            if equalFlagGeom && equalFlag0 && equalFlag1 && equalFlag2 && equalFlag3 && equalFlag4 && equalFlag5
                mainStruct.paramLogStruct(k1,1).singleCrackNo(1) = k2;
            elseif k2 == length(singleStruct.paramLogStruct);
                error(['Cannot find a crack with the same dimensions as Crack A in singleStruct. Crack pair no. ',num2str(k1)])
            end
            k2 = k2+1;
        end
        
        %For Crack B
        k2 = 1;
        while mainStruct.paramLogStruct(k1,1).singleCrackNo(2) == 0
            %Dynamically unpack singleStruct.paramLogStruct(k2,1).naturalParams
            s = struct;
            fieldNames = fieldnames(singleStruct.paramLogStruct(k2,1).naturalParams);
            for j1 = 1:length(fieldNames)
                if isstruct(eval(['singleStruct.paramLogStruct(k2,1).naturalParams.',fieldNames{j1},';']))
                    subFieldNames = fieldnames(eval(['singleStruct.paramLogStruct(k2,1).naturalParams.',fieldNames{j1},';']));
                    for j2 = 1:length(subFieldNames)
                        if isfield(s,subFieldNames{j2})
                            error(['Error unpacking singleStruct.paramLogStruct(k2,1).naturalParams. Duplicate field name: ',subFieldNames{j2}]);
                        else
                            eval(['s.',subFieldNames{j2},' = singleStruct.paramLogStruct(k2,1).naturalParams.',fieldNames{j1},'.',subFieldNames{j2},';']);
                        end
                    end
                else
                    if isfield(s,fieldNames{j1})
                        error(['Error unpacking singleStruct.paramLogStruct(k2,1).naturalParams. Duplicate field name: ',fieldNames{j1}]);
                    else
                        eval(['s.',fieldNames{j1},' = singleStruct.paramLogStruct(k2,1).naturalParams.',fieldNames{j1},';']);
                    end
                end
            end
            if twoSingleFoldersFlag   %If there is a separate set of models for Crack B, load singleStruct2 for identical cracks.
                %Dynamically unpack singleStruct2.paramLogStruct(k2,1).naturalParams
                s2 = struct;
                fieldNames = fieldnames(singleStruct2.paramLogStruct(k2,1).naturalParams);
                for j1 = 1:length(fieldNames)
                    if isstruct(eval(['singleStruct2.paramLogStruct(k2,1).naturalParams.',fieldNames{j1},';']))
                        subFieldNames = fieldnames(eval(['singleStruct2.paramLogStruct(k2,1).naturalParams.',fieldNames{j1},';']));
                        for j2 = 1:length(subFieldNames)
                            if isfield(s2,subFieldNames{j2})
                                error(['Error unpacking singleStruct2.paramLogStruct(k2,1).naturalParams. Duplicate field name: ',subFieldNames{j2}]);
                            else
                                eval(['s2.',subFieldNames{j2},' = singleStruct2.paramLogStruct(k2,1).naturalParams.',fieldNames{j1},'.',subFieldNames{j2},';']);
                            end
                        end
                    else
                        if isfield(s2,fieldNames{j1})
                            error(['Error unpacking singleStruct2.paramLogStruct(k2,1).naturalParams. Duplicate field name: ',fieldNames{j1}]);
                        else
                            eval(['s2.',fieldNames{j1},' = singleStruct2.paramLogStruct(k2,1).naturalParams.',fieldNames{j1},';']);
                        end
                    end
                end                
                
                %Check for model equivalence
                equalFlagGeom = strcmpi(m.geometryType,s2.geometryType);        %General Geometry Type
                if strcmpi(m.geometryType,'pipe')
                    if m.oppositeFlag && ~m.subsurfaceBFlag
                        equalFlagGeom = equalFlagGeom && ~strcmpi(m.pipeCrackA,s2.pipeCrackA);	%If Crack B is opposite Crack A, the single crack model should be for the opposite pipe surface.
                    else
                        equalFlagGeom = equalFlagGeom && strcmpi(m.pipeCrackA,s2.pipeCrackA);	%...otherwise Crack B should be on same pipe surface.
                    end
                    equalFlag0 = isequal(m.ri,s2.ri) && isequal(m.b,s2.b);          %Pipe radius & wall thickness
                else
                    equalFlag0 = isequal(m.b,s2.b);                                 %Plate thickness
                end
                equalFlag1 = abs(mainStruct.paramLogStruct(k1,1).modelParams.crackSizes.aB1-...
                    singleStruct2.paramLogStruct(k2,1).modelParams.crackSizes.aA1)...
                    < 1e-3*singleStruct2.paramLogStruct(k2,1).modelParams.crackSizes.aA1;   %Crack size in y dir. - note use of modelParams
                %Subsurface crack aspect ratio - in mainStruct aB1OveraB2
                %can be NaN indicating that it is equal to aA1OveraA2.
                if isnan(m.aB1OveraB2)
                    equalFlag2 = abs(m.aA1OveraA2-s2.aA1OveraA2) < 1e-3*s2.aA1OveraA2;
                else
                    equalFlag2 = abs(m.aB1OveraB2-s2.aA1OveraA2) < 1e-3*s2.aA1OveraA2;
                end
                equalFlag3 = isequal(m.nu,s2.nu);                             %Poisson's ratio
                equalFlag4 = isequal(m.subsurfaceBFlag,s2.subsurfaceAFlag);   %Subsurface or not
                %Subsurface crack depth
                if m.subsurfaceBFlag
                    equalFlag5 = abs(mainStruct.paramLogStruct(k1,1).modelParams.crackPositions.yB-...
                        singleStruct2.paramLogStruct(k2,1).modelParams.crackPositions.yA)...
                        < 1e-3*singleStruct2.paramLogStruct(k2,1).modelParams.crackPositions.yA;    %Subsurface crack depth - note use of modelParams
                else
                    equalFlag5 = true;
                end
                %Loading conditions
                if strcmpi(m.loadingType,'crackFacePressure')
                    equalFlag6 = strcmpi(m.CFPstring,s2.CFPstring);
                elseif any(strcmpi(m.loadingType,{'pipePressure','pipePressureInclCrack'}))
                    equalFlag6 = strcmpi(m.loadingType,s2.loadingType) && isequal(m.pressureMagnitude,s2.pressureMagnitude);
                else
                    equalFlag6 = strcmpi(m.loadingType,s2.loadingType) && isequal(m.loadMagnitude,s2.loadMagnitude);
                end
                %Loading conditions - biaxial direction
                equalFlag7 = isequal(m.biaxialFlag,s2.biaxialFlag) && (isequal(m.biaxialStress,s2.biaxialStress) || (isnan(m.biaxialStress) && isnan(s2.biaxialStress)));  %Loading conditions - biaxial direction
                %Check that all parameters are (simultaneously) correct
                if equalFlagGeom && equalFlag0 && equalFlag1 && equalFlag2 && equalFlag3 && equalFlag4 && equalFlag5 && equalFlag6 && equalFlag7
                    mainStruct.paramLogStruct(k1,1).singleCrackNo(2) = k2;
                elseif k2 == length(singleStruct2.paramLogStruct);
                    error(['Cannot find a crack with the same dimensions as Crack B in singleStruct2. Crack pair no. ',num2str(k1)])
                end
            else                %(i.e. If there is only one folder containing single defect results)
                %Check for model equivalence
                equalFlagGeom = strcmpi(m.geometryType,s.geometryType);                 %General Geometry Type
                if strcmpi(m.geometryType,'pipe')
                    equalFlagGeom = equalFlagGeom && strcmpi(m.pipeCrackA,s.pipeCrackA);
                    equalFlag0 = isequal(m.ri,s.ri) && isequal(m.b,s.b);                    %Pipe radius & wall thickness
                else
                    equalFlag0 = isequal(m.b,s.b);                                          %Plate thickness
                end
                equalFlag1 = abs(mainStruct.paramLogStruct(k1,1).modelParams.crackSizes.aB1-...
                    singleStruct.paramLogStruct(k2,1).modelParams.crackSizes.aA1)...
                    < 1e-3*singleStruct.paramLogStruct(k2,1).modelParams.crackSizes.aA1;    %Crack depth - note use of modelParams
                equalFlag2 = abs(m.aB1OveraB2-s.aA1OveraA2) < 1e-3*s.aA1OveraA2;            %Crack aspect ratio
                equalFlag3 = isequal(m.nu,s.nu);              %Poisson's ratio
                %Loading conditions
                if strcmpi(m.loadingType,'crackFacePressure')
                    equalFlag4 = strcmpi(m.CFPstring,s.CFPstring);
                elseif any(strcmpi(m.loadingType,{'pipePressure','pipePressureInclCrack'}))
                    equalFlag4 = strcmpi(m.loadingType,s.loadingType) && isequal(m.pressureMagnitude,s.pressureMagnitude);
                else
                    equalFlag4 = strcmpi(m.loadingType,s.loadingType) && isequal(m.loadMagnitude,s.loadMagnitude);
                end
                %Loading conditions - biaxial direction
                equalFlag5 = isequal(m.biaxialFlag,s.biaxialFlag) && (isequal(m.biaxialStress,s.biaxialStress) || (isnan(m.biaxialStress) && isnan(s.biaxialStress)));  %Loading conditions - biaxial direction
                %Check that all parameters are (simultaneously) correct
                if equalFlagGeom && equalFlag0 && equalFlag1 && equalFlag2 && equalFlag3 && equalFlag4 && equalFlag5
                    mainStruct.paramLogStruct(k1,1).singleCrackNo(2) = k2;
                elseif k2 == length(singleStruct.paramLogStruct);
                    error(['Cannot find a crack with the same dimensions as Crack B in singleStruct. Crack pair no. ',num2str(k1)])
                end
            end
            k2 = k2+1;
        end
        
        %Assign this to an additional variable to make the code more readable.
        singleCrackNo = mainStruct.paramLogStruct(k1,1).singleCrackNo;
        
        %Calculate interaction factors and add them to mainStruct.outputArray
        if mainStruct.outputArray{1,k1}.datReadSuccessFlag &&...
                mainStruct.paramLogStruct(k1).jobExitStatus &&...
                ~mainStruct.paramLogStruct(k1).contDepWarnFlagLevel2 &&...
                ~isempty(mainStruct.outputArray{1,k1}.k{2,1}) &&...
                ~isempty(mainStruct.outputArray{1,k1}.k{2,1}) %If the main result is thought to be okay
            
            %Calculate local KI interaction factors
            for k3 = [1,2]  %Note you now need this additional loop for calculating local interaction factors
                %First create singleResultValidFlag, which indicates
                %whether we think the result for a single crack for this
                %case is okay.
                if twoSingleFoldersFlag
                    if k3==1
                        singleResultValidFlag = singleStruct.outputArray{1,singleCrackNo(k3)}.datReadSuccessFlag &&...
                            singleStruct.paramLogStruct(singleCrackNo(k3)).jobExitStatus &&...
                            ~singleStruct.paramLogStruct(singleCrackNo(k3)).contDepWarnFlagLevel2 &&...
                            ~isempty(singleStruct.outputArray{1,singleCrackNo(k3)}.k{2,1});
                    else
                        singleResultValidFlag = singleStruct2.outputArray{1,singleCrackNo(k3)}.datReadSuccessFlag &&...
                            singleStruct2.paramLogStruct(singleCrackNo(k3)).jobExitStatus &&...
                            ~singleStruct2.paramLogStruct(singleCrackNo(k3)).contDepWarnFlagLevel2 &&...
                            ~isempty(singleStruct2.outputArray{1,singleCrackNo(k3)}.k{2,1});
                    end
                else
                    singleResultValidFlag = singleStruct.outputArray{1,singleCrackNo(k3)}.datReadSuccessFlag &&...
                        singleStruct.paramLogStruct(singleCrackNo(k3)).jobExitStatus &&...
                        ~singleStruct.paramLogStruct(singleCrackNo(k3)).contDepWarnFlagLevel2 &&...
                        ~isempty(singleStruct.outputArray{1,singleCrackNo(k3)}.k{2,1});
                end
                
                %Now, if the single-crack result is thought to be okay,
                %calculate the interaction factors.
                if singleResultValidFlag
                    try
                        if twoSingleFoldersFlag     %Crack B is a subsurface crack.
                            if k3 == 1
                                mainStruct.outputArray{1,k1}.interactionFactor{1,k3} = mainStruct.outputArray{1,k1}.k{2,k3}./singleStruct.outputArray{1,singleCrackNo(k3)}.k{2,1};
                            else
                                %NOTE - The order of crack tip points on the single crack is reversed so that it matches with the order of points on the
                                %embedded crack in the proximate crack model.
                                mainStruct.outputArray{1,k1}.interactionFactor{1,k3} = mainStruct.outputArray{1,k1}.k{2,k3}./fliplr(singleStruct2.outputArray{1,singleCrackNo(k3)}.k{2,1});
                            end
                        else    %Crack B is a surface crack.
                            mainStruct.outputArray{1,k1}.interactionFactor{1,k3} = mainStruct.outputArray{1,k1}.k{2,k3}./singleStruct.outputArray{1,singleCrackNo(k3)}.k{2,1};
                        end
                        mainStruct.outputArray{1,k1}.singleIntFactorsCalced(k3) = true;
                    catch
                        warning(['Problem calculating local interaction factors for crack pair ',num2str(k1),', Crack ',num2str(k3),'.']);
                        mainStruct.outputArray{1,k1}.singleIntFactorsCalced(k3) = false;
                        mainStruct.outputArray{1,k1}.intFactorsCalced = false;
                    end
                else
                    warning(['Interaction factors not calculated for model no. ',num2str(k1),' crack ',num2str(k3),' no valid result for single crack.'])
                    mainStruct.outputArray{1,k1}.interactionFactor{1,k3} = [];
                    mainStruct.outputArray{1,k1}.singleIntFactorsCalced(k3) = false;
                    mainStruct.outputArray{1,k1}.intFactorsCalced = false;
                end
            end
            
            %Calculate global KI interaction factors
            if ~any(~mainStruct.outputArray{1,k1}.singleIntFactorsCalced)  %This bit is only executed if you have local interaction factors for both cracks
                %Global interaction factors
                if twoSingleFoldersFlag   %Crack B is a subsurface crack.
                    %Calculate global interaction factors
                    mainStruct.outputArray{1,k1}.interactionFactorGlobal{1,1} = ...
                        max(mainStruct.outputArray{1,k1}.k{2,1},[],2)./max(singleStruct.outputArray{1,singleCrackNo(1)}.k{2,1},[],2);  %Using maximum KI over A only
                    
                    mainStruct.outputArray{1,k1}.interactionFactorGlobal{1,2} = ...
                        max(mainStruct.outputArray{1,k1}.k{2,2},[],2)./max(singleStruct2.outputArray{1,singleCrackNo(2)}.k{2,1},[],2);  %Using maximum KI over B only
                    
                    mainStruct.outputArray{1,k1}.interactionFactorGlobal{1,3} = ...
                        max([mainStruct.outputArray{1,k1}.k{2,1},mainStruct.outputArray{1,k1}.k{2,2}],[],2)./...
                        max([singleStruct.outputArray{1,singleCrackNo(1)}.k{2,1},singleStruct2.outputArray{1,singleCrackNo(2)}.k{2,1}],[],2);  %Using maximum KI over both A and B
                    
                    %Locations for the global interaction factors
                    [~,mainStruct.outputArray{1,k1}.interactionFactorGlobal{2,1}] = max(mainStruct.outputArray{1,k1}.k{2,1},[],2);  %A
                    [~,mainStruct.outputArray{1,k1}.interactionFactorGlobal{3,1}] = max(singleStruct.outputArray{1,singleCrackNo(1)}.k{2,1},[],2);
                    %
                    [~,mainStruct.outputArray{1,k1}.interactionFactorGlobal{2,2}] = max(mainStruct.outputArray{1,k1}.k{2,2},[],2);  %B
                    [~,mainStruct.outputArray{1,k1}.interactionFactorGlobal{3,2}] = max(singleStruct2.outputArray{1,singleCrackNo(2)}.k{2,1},[],2);
                else   %Crack B is a subsurface crack.
                    %Calculate global interaction factors
                    mainStruct.outputArray{1,k1}.interactionFactorGlobal{1,1} = ...
                        max(mainStruct.outputArray{1,k1}.k{2,1},[],2)./max(singleStruct.outputArray{1,singleCrackNo(1)}.k{2,1},[],2);  %Using maximum KI over A only
                    
                    mainStruct.outputArray{1,k1}.interactionFactorGlobal{1,2} = ...
                        max(mainStruct.outputArray{1,k1}.k{2,2},[],2)./max(singleStruct.outputArray{1,singleCrackNo(2)}.k{2,1},[],2);  %Using maximum KI over B only
                    
                    mainStruct.outputArray{1,k1}.interactionFactorGlobal{1,3} = ...
                        max([mainStruct.outputArray{1,k1}.k{2,1},mainStruct.outputArray{1,k1}.k{2,2}],[],2)./...
                        max([singleStruct.outputArray{1,singleCrackNo(1)}.k{2,1},singleStruct.outputArray{1,singleCrackNo(2)}.k{2,1}],[],2);  %Using maximum KI over both A and B
                    
                    %Locations for the global interaction factors
                    [~,mainStruct.outputArray{1,k1}.interactionFactorGlobal{2,1}] = max(mainStruct.outputArray{1,k1}.k{2,1},[],2);  %A
                    [~,mainStruct.outputArray{1,k1}.interactionFactorGlobal{3,1}] = max(singleStruct.outputArray{1,singleCrackNo(1)}.k{2,1},[],2);
                    %
                    [~,mainStruct.outputArray{1,k1}.interactionFactorGlobal{2,2}] = max(mainStruct.outputArray{1,k1}.k{2,2},[],2);  %B
                    [~,mainStruct.outputArray{1,k1}.interactionFactorGlobal{3,2}] = max(singleStruct.outputArray{1,singleCrackNo(2)}.k{2,1},[],2);
                end
                
                %Save a logical to indicate whether IFs were calculated for this model
                mainStruct.outputArray{1,k1}.intFactorsCalced = true;
            end
        else
            warning(['Interaction factors not calculated for model no. ',num2str(k1),' no valid result for the crack pair.'])
            mainStruct.outputArray{1,k1}.interactionFactor{1,1} = [];
            mainStruct.outputArray{1,k1}.singleIntFactorsCalced = [false,false];
            mainStruct.outputArray{1,k1}.intFactorsCalced = false;
        end
    end
else
    for k1 = 1:length(mainStruct.outputArray)
        %Display the model number for every 100th model
        if mod(k1,100) == 0
            disp(['Current model no.: ',num2str(k1)])
        end
        
        %Dynamically unpack mainStruct.paramLogStruct(k1,1).naturalParams
        m = struct;
        fieldNames = fieldnames(mainStruct.paramLogStruct(k1,1).naturalParams);
        for j1 = 1:length(fieldNames)
            if isstruct(eval(['mainStruct.paramLogStruct(k1,1).naturalParams.',fieldNames{j1},';']))
                subFieldNames = fieldnames(eval(['mainStruct.paramLogStruct(k1,1).naturalParams.',fieldNames{j1},';']));
                for j2 = 1:length(subFieldNames)
                    if isfield(m,subFieldNames{j2})
                        error(['Error unpacking mainStruct.paramLogStruct(k1,1).naturalParams. Duplicate field name: ',subFieldNames{j2}]);
                    else
                        eval(['m.',subFieldNames{j2},' = mainStruct.paramLogStruct(k1,1).naturalParams.',fieldNames{j1},'.',subFieldNames{j2},';']);
                    end
                end
            else
                if isfield(m,fieldNames{j1})
                    error(['Error unpacking mainStruct.paramLogStruct(k1,1).naturalParams. Duplicate field name: ',fieldNames{j1}]);
                else
                    eval(['m.',fieldNames{j1},' = mainStruct.paramLogStruct(k1,1).naturalParams.',fieldNames{j1},';']);
                end
            end
        end
        
        %Work out which model for a single crack to use from the naturalParams
        mainStruct.paramLogStruct(k1,1).singleCrackNo = [];
        k2 = 1;
        while isempty(mainStruct.paramLogStruct(k1,1).singleCrackNo)
            %Dynamically unpack singleStruct.paramLogStruct(k2,1).naturalParams
            s = struct;
            fieldNames = fieldnames(singleStruct.paramLogStruct(k2,1).naturalParams);
            for j1 = 1:length(fieldNames)
                if isstruct(eval(['singleStruct.paramLogStruct(k2,1).naturalParams.',fieldNames{j1},';']))
                    subFieldNames = fieldnames(eval(['singleStruct.paramLogStruct(k2,1).naturalParams.',fieldNames{j1},';']));
                    for j2 = 1:length(subFieldNames)
                        if isfield(s,subFieldNames{j2})
                            error(['Error unpacking singleStruct.paramLogStruct(k2,1).naturalParams. Duplicate field name: ',subFieldNames{j2}]);
                        else
                            eval(['s.',subFieldNames{j2},' = singleStruct.paramLogStruct(k2,1).naturalParams.',fieldNames{j1},'.',subFieldNames{j2},';']);
                        end
                    end
                else
                    if isfield(s,subFieldNames{j1})
                        error(['Error unpacking singleStruct.paramLogStruct(k2,1).naturalParams. Duplicate field name: ',fieldNames{j1}]);
                    else
                        eval(['s.',fieldNames{j1},' = singleStruct.paramLogStruct(k2,1).naturalParams.',fieldNames{j1},';']);
                    end
                end
            end
            
            %Check for model equivalence
            equalFlag1 = isequal(m.aA1Overb,s.aA1Overb);
            equalFlag2 = abs(m.aA1OveraA2-s.aA1OveraA2) < 1e-3*s.aA1OveraA2;
            equalFlag3 = isequal(m.aB1OveraA1,s.aB1OveraA1);
            equalFlag4 = abs(m.aB1OveraB2-s.aB1OveraB2) < 1e-3*s.aB1OveraB2;
            if isnan(m.aB1OveraB2) && isnan(s.aB1OveraB2)
                equalFlag4 = true;
            end
            equalFlag5 = isequal(m.nu,s.nu);                    %Poisson's ratio
            %Loading conditions
            if strcmpi(m.loadingType,'crackFacePressure')
                equalFlag6 = strcmpi(m.CFPstring,s.CFPstring);
            else
                equalFlag6 = strcmpi(m.loadingType,s.loadingType) && isequal(m.loadMagnitude,s.loadMagnitude);
            end
            %Check that all parameters are (simultaneously) correct
            if equalFlag1 && equalFlag2 && equalFlag3 && equalFlag4 && equalFlag5 && equalFlag6
                mainStruct.paramLogStruct(k1,1).singleCrackNo = k2;
            end
            k2 = k2+1;
        end
        
        %Assign this to a additional variable to make code more readable
        singleCrackNo = mainStruct.paramLogStruct(k1,1).singleCrackNo;
        
        %Calculate interaction factors and add them to mainStruct.outputArray
        if mainStruct.outputArray{1,k1}.datReadSuccessFlag &&...
                mainStruct.paramLogStruct(k1).jobExitStatus &&...
                ~mainStruct.paramLogStruct(k1).contDepWarnFlagLevel2 &&...
                ~isempty(mainStruct.outputArray{1,k1}.k{2,1}) &&...
                ~isempty(mainStruct.outputArray{1,k1}.k{2,1}) %If the main result is thought to be okay
            if singleStruct.outputArray{1,singleCrackNo}.datReadSuccessFlag &&...
                    singleStruct.paramLogStruct(singleCrackNo).jobExitStatus &&...
                    ~singleStruct.paramLogStruct(singleCrackNo).contDepWarnFlagLevel2 &&...
                    ~isempty(singleStruct.outputArray{1,singleCrackNo}.k{2,1}) &&...
                    ~isempty(singleStruct.outputArray{1,singleCrackNo}.k{2,2})   %If the single-crack result is thought to be okay
                %Calculate KI interaction factors
                %Local
                mainStruct.outputArray{1,k1}.interactionFactor{1,1} = mainStruct.outputArray{1,k1}.k{2,1}./singleStruct.outputArray{1,singleCrackNo}.k{2,1};    %A
                mainStruct.outputArray{1,k1}.interactionFactor{1,2} = mainStruct.outputArray{1,k1}.k{2,2}./singleStruct.outputArray{1,singleCrackNo}.k{2,2};    %B
                
                %'Global'
                mainStruct.outputArray{1,k1}.interactionFactorGlobal{1,1} = ...
                    max(mainStruct.outputArray{1,k1}.k{2,1},[],2)./max(singleStruct.outputArray{1,singleCrackNo}.k{2,1},[],2);  %Using maximum KI over A only
                mainStruct.outputArray{1,k1}.interactionFactorGlobal{1,2} = ...
                    max(mainStruct.outputArray{1,k1}.k{2,2},[],2)./max(singleStruct.outputArray{1,singleCrackNo}.k{2,2},[],2);  %Using maximum KI over B only
                mainStruct.outputArray{1,k1}.interactionFactorGlobal{1,3} = ...
                    max([mainStruct.outputArray{1,k1}.k{2,1},mainStruct.outputArray{1,k1}.k{2,2}],[],2)./...
                    max([singleStruct.outputArray{1,singleCrackNo}.k{2,1},singleStruct.outputArray{1,singleCrackNo}.k{2,2}],[],2);  %Using maximum KI over both A and B
                
                %Locations for the global interaction factors
                [~,mainStruct.outputArray{1,k1}.interactionFactorGlobal{2,1}] = max(mainStruct.outputArray{1,k1}.k{2,1},[],2);  %A
                [~,mainStruct.outputArray{1,k1}.interactionFactorGlobal{3,1}] = max(singleStruct.outputArray{1,singleCrackNo}.k{2,1},[],2);
                
                [~,mainStruct.outputArray{1,k1}.interactionFactorGlobal{2,2}] = max(mainStruct.outputArray{1,k1}.k{2,2},[],2);  %B
                [~,mainStruct.outputArray{1,k1}.interactionFactorGlobal{3,2}] = max(singleStruct.outputArray{1,singleCrackNo}.k{2,2},[],2);
                
                %Save a logical to indicate whether IFs were calculated for this model
                mainStruct.outputArray{1,k1}.intFactorsCalced = true;
            else
                warning(['Interaction factors not calculated for model no. ',num2str(k1),' no valid result for single crack.'])
                mainStruct.outputArray{1,k1}.interactionFactor{1,1} = [];
                mainStruct.outputArray{1,k1}.intFactorsCalced = false;
            end
        else
            warning(['Interaction factors not calculated for model no. ',num2str(k1),' no valid result for the crack pair.'])
            mainStruct.outputArray{1,k1}.interactionFactor{1,1} = [];
            mainStruct.outputArray{1,k1}.intFactorsCalced = false;
        end
    end
end
disp('DONE!')

%% Calculate how many of the models the interaction factor was calculated for and display
totalIntFactorsCalced = 0;
for k1 = 1:length(mainStruct.outputArray)
    if mainStruct.outputArray{1,k1}.intFactorsCalced
        totalIntFactorsCalced = totalIntFactorsCalced + 1;
    end
end
disp(['Interaction factors calculated for ',num2str(totalIntFactorsCalced),' out of ',num2str(length(mainStruct.outputArray)),' models.']);

%% Save results file
disp('Saving results file...')
save('results.mat');

%% Write logfile
disp('Writing logfile...')
if exist('singleStruct2','var')
    int_defects_write_logfile_2( mainStruct, singleStruct, singleStruct2);
else
    int_defects_write_logfile_2( mainStruct, singleStruct);
end

disp('DONE!')

end