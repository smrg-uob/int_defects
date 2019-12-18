function [ outputArray, varargout ] = read_dat_contour_integral( inputStruct, contourType, varargin )
%read_dat_contour_integral.m
%Harry Coules 2015
%
%DESCRIPTION
%This function is used to read contour integral output data from an Abaqus
%.dat file.
%
%INPUT ARGUMENTS
%   inputStruct - Structure containing the name of the .dat file
%   (inputStruct.filename). Alternatively, if the first input argument is a
%   string, this is taken to be the filename.
%   contourType - A string specifying the type of contour integral output.
%   It can either be 'k' (stress intensity factors), 't' (T-stress) or 'j'
%   (J-integral).
%*OPTIONAL INPUT ARGUMENT*
%   stepsToRead - Model steps for which output is to be read, as a vector.
%    If this is omitted, output is read for all steps.
%
%OUTPUT ARGUMENTS
%   outputArray - cell array contining either stress intensity factors,
%       J-integral values, or T-stress values, depending on the type of contour
%       integral results present (defined by contourType).
%*OPTIONAL OUTPUT ARGUMENT*
%   stepIncTimes - Array containing information on the steps and increments
%   found in the .dat file. It is an n-by-6 array with the columns: step no.,
%   inc. no., time increment completed, fraction of step completed, step
%   time completed, total time completed.
%
%NOTES
% - In this function, the model steps to read are specified as an
%   additional optional input variable. Note that this is different to the
%   method used in read_dat.m, where the steps are specifed as a field in
%   inputStruct.

%% Process input
%Check the first input argument...
%It can be a struct containing the field .filename OR it can be a string
%which alone defines the filename
if isstruct(inputStruct)  %If it's a struct
    %If filename is empty in input structure, give an error.
    if ~isfield(inputStruct,'filename')
        error('No filename defined for the .dat file to be read.');
    else
        filename = inputStruct.filename;
    end
elseif ischar(inputStruct)    %If it's a string
    filename = inputStruct;
end

%Check the second input argument...
if ~strcmpi(contourType,'j') && ~strcmpi(contourType,'k') && ~strcmpi(contourType,'t') && ~strcmpi(contourType,'none')
    error('Input argument contourType must be either ''j'', ''k'', ''t'', or ''none''.')
end

%Optional input argument: stepsToRead
if ~isempty(varargin)
    if length(varargin) == 1
        stepsToRead = varargin{1};
    else
        error('Too many input arguments.');
    end
else
    stepsToRead = [];
end

%% Read in the lines of the specified file
fid = fopen(filename,'r','n','UTF-8');
C = textscan(fid, '%s', Inf, 'Delimiter','\n'); C = C{1};
fclose(fid);

%% Determine which lines indicate the start of a step output sequence
stepHeadersLogical = [];
stepHeadersIndex = [];
stepNos = [];

incHeadersLogical = [];
incHeadersIndex = [];
stepIncTimes = [];

for k1 = 1:length(C)
    %Logical indicating the lines in which the step headers occur
    stepHeadersLogical(k1) = ~isempty(sscanf(C{k1},'S T E P       %d     S T A T I C   A N A L Y S I S'))...
        && isempty(strfind(C{k1},'INCREMENT'));
    
    if stepHeadersLogical(k1)
        stepHeadersIndex = [stepHeadersIndex;k1];   %Indices of the lines in which the step headers occur
        stepNos = [stepNos;sscanf(C{k1},'S T E P       %d     S T A T I C   A N A L Y S I S')]; %Step numbers
    end
    
    %Logical indicating the lines in which the increment headers occur
    incHeadersLogical(k1) = ~isempty(strfind(C{k1},'INCREMENT')) && ~isempty(strfind(C{k1},'SUMMARY'));
    
    if incHeadersLogical(k1)
        incHeadersIndex = [incHeadersIndex;k1]; %Indices of the lines in which the increment headers occur
        stepIncTimes = [stepIncTimes;...
            stepNos(end),...
            sscanf(C{k1},'INCREMENT     %d SUMMARY'),...
            sscanf(C{k1+3},'TIME INCREMENT COMPLETED %f , FRACTION OF STEP COMPLETED %f')',...
            sscanf(C{k1+4},'STEP TIME COMPLETED %f , TOTAL TIME COMPLETED %f')']; %Array with: Step no, inc no, time inc completed, fraction of step completed, step time completed, total time completed.
    end
end

%% Determine which lines contain the start of a contour integral output sequence
kHeadersLogical = strcmpi('K   F A C T O R       E S T I M A T E S',C);
jHeadersLogical = strcmpi('J - I N T E G R A L   E S T I M A T E S',C);
tHeadersLogical = strcmpi('T - S T R E S S   E S T I M A T E S',C);

kHeadersIndex = find(kHeadersLogical);
jHeadersIndex = find(jHeadersLogical);
tHeadersIndex = find(tHeadersLogical);

%Determine which step each contour integral output sequence belongs to
kStep = [];
for k1 = 1:length(kHeadersIndex)
    afterStepHeaderLogical = stepHeadersIndex < kHeadersIndex(k1);
    [~,stepIndex] = max(stepHeadersIndex(afterStepHeaderLogical));
    kStep(k1) = stepNos(stepIndex);
end

jStep = [];
for k1 = 1:length(jHeadersIndex)
    afterStepHeaderLogical = stepHeadersIndex < jHeadersIndex(k1);
    [~,stepIndex] = max(stepHeadersIndex(afterStepHeaderLogical));
    jStep(k1) = stepNos(stepIndex);
end

tStep = [];
for k1 = 1:length(tHeadersIndex)
    afterStepHeaderLogical = stepHeadersIndex < tHeadersIndex(k1);
    [~,stepIndex] = max(stepHeadersIndex(afterStepHeaderLogical));
    tStep(k1) = stepNos(stepIndex);
end

%% Extract the data from the contour integral output sequence(s)
if strcmpi(contourType,'k')
    %% Extract data for k contour integral from output sequence(s)
    for n1 = 1:length(kHeadersIndex)
        %Create cell arrays
        kContoursCell{n1} = [];
        k1Cell{n1} = [];
        k2Cell{n1} = [];
        k3Cell{n1} = [];
        kMTSDirCell{n1} = [];
        jFromKCell{n1} = [];
        n2 = kHeadersIndex(n1);
        endFlag = false;
        while ~endFlag
            n2 = n2+1;
            if strcmpi(char(C(n2)),'LABELS REFERENCED IN THE ABOVE TABLE');
                outputFlag = 'LABELS';   %Data labels
                endFlag = true;
            elseif strcmpi(char(C(n2)),'THE ANALYSIS HAS BEEN COMPLETED');
                outputFlag = 'END'; %End of the .dat file
                endFlag = true;
            else
                %Special regexp for K factor output
                if strfind(char(C(n2)), 'K1:')
                    if strfind(char(C(n2)), 'XFEM') %This is for the case of SIF output from an XFEM analysis
                        k1Line = sscanf(char(C(n2)),' XFEM_%u       K1: %f %f %f %f %f');
                        k1Line(1) = [];
                    else
                        k1Line = sscanf(char(C(n2)),' K1: %f %f %f %f %f'); %SIF output from a normal fracture analysis 
                        if isempty(k1Line)
                            k1Line = sscanf(char(C(n2)),' -%u- K1: %f %f %f %f %f');
                            k1Line(1) = [];
                        end
                    end
                    k1Cell{n1} = [k1Cell{n1};k1Line];
                elseif strfind(char(C(n2)), 'K2:')
                    k2Line = sscanf(char(C(n2)),' K2: %f %f %f %f %f');
                    k2Cell{n1} = [k2Cell{n1};k2Line];
                elseif strfind(char(C(n2)), 'K3:')
                    k3Line = sscanf(char(C(n2)),' K3: %f %f %f %f %f');
                    k3Cell{n1} = [k3Cell{n1};k3Line];
                elseif strfind(char(C(n2)), 'MTS   DIRECTION (DEG):')
                    kMTSDirLine = sscanf(char(C(n2)),' MTS   DIRECTION (DEG): %f %f %f %f %f');
                    kMTSDirCell{n1} = [kMTSDirCell{n1};kMTSDirLine];
                elseif strfind(char(C(n2)), 'J from Ks:')
                    jFromKLine = sscanf(char(C(n2)),' J from Ks: %f %f %f %f %f');
                    jFromKCell{n1} = [jFromKCell{n1};jFromKLine];
                else kContoursLine = sscanf(char(C(n2)),'%f %f %f %f %f');
                    kContoursCell{n1} = [kContoursCell{n1};kContoursLine];
                end
            end
        end
        k1Cell{n1} = reshape(k1Cell{n1},length(kContoursCell{n1}),[]);    %Reshape arrays. Contour no down colums, point on crack along rows.
        k2Cell{n1} = reshape(k2Cell{n1},length(kContoursCell{n1}),[]);
        k3Cell{n1} = reshape(k3Cell{n1},length(kContoursCell{n1}),[]);
        kMTSDirCell{n1} = reshape(kMTSDirCell{n1},length(kContoursCell{n1}),[]);
        jFromKCell{n1} = reshape(jFromKCell{n1},length(kContoursCell{n1}),[]);
    end
elseif strcmpi(contourType,'j')
    %% Extract data for j contour integral from output sequence(s)
    for n1 = 1:length(jHeadersIndex)
        jInt1Cell{n1} = [];
        jInt2Cell{n1} = [];
        n2 = jHeadersIndex(n1);
        endFlag = false;
        dataSwitch = false;
        while ~endFlag
            n2 = n2+1;
            if strcmpi(char(C(n2)),'LABELS REFERENCED IN THE ABOVE TABLE');
                outputFlag = 'LABELS';   %Data labels
                endFlag = true;
            else
                %Special regexp for J-integral output format
                if ~isempty(strfind(char(C(n2)), '- '))
                    jLine = sscanf(char(C(n2)),' -%u- %f %f %f %f %f')';
                    if length(jLine) > 1
                        jLine = jLine(2:end);
                        dataSwitch = true; %If there's a -%u- in the line, and at least one other number, we've moved on to the J-Integral data lines
                    end
                else
                    jLine = sscanf(char(C(n2)),'%f %f %f %f %f')';
                end
                if ~isempty(jLine)
                    if ~exist('prevJLineLength','var')
                        prevJLineLength = Inf;
                        dataSwitch = false;
                    end
                    if ~dataSwitch;
                        jInt1Cell{n1} = [jInt1Cell{n1};jLine'];
                    else
                        jInt2Cell{n1} = [jInt2Cell{n1};jLine'];
                    end
                    prevJLineLength = length(jLine);
                end
            end
        end
        jInt2Cell{n1} = reshape(jInt2Cell{n1},length(jInt1Cell{n1}),[]);    %Reshape array. Contour no down colums, point on crack along rows.
    end
elseif strcmpi(contourType,'t')
    %% Extract data for t contour integral from output sequence(s)
    for n1 = 1:length(tHeadersIndex)
        tStress1Cell{n1} = [];
        tStress2Cell{n1} = [];
        n2 = tHeadersIndex(n1);
        endFlag = false;
        dataSwitch = false;
        while ~endFlag
            n2 = n2+1;
            if strcmpi(char(C(n2)),'LABELS REFERENCED IN THE ABOVE TABLE');
                outputFlag = 'LABELS';   %Data labels
                endFlag = true;
            else
                %Special regexp for T-stress output format (note it's the same as for the J integral)
                if ~isempty(strfind(char(C(n2)), '- '))
                    tLine = sscanf(char(C(n2)),' -%u- %f %f %f %f %f')';
                    if length(tLine) > 1
                        tLine = tLine(2:end);
                        dataSwitch = true; %If there's a -%u- in the line, and at least one other number, we've moved on to the T-stress data lines
                    end
                else
                    tLine = sscanf(char(C(n2)),'%f %f %f %f %f')';
                end
                if ~isempty(tLine)
                    if ~exist('prevTLineLength','var')
                        prevTLineLength = Inf;
                        dataSwitch = false;
                    end
                    if ~dataSwitch
                        tStress1Cell{n1} = [tStress1Cell{n1};tLine'];
                    else
                        tStress2Cell{n1} = [tStress2Cell{n1};tLine'];
                    end
                    prevTLineLength = length(tLine);
                end
            end
        end
        tStress2Cell{n1} = reshape(tStress2Cell{n1},length(tStress1Cell{n1}),[]);    %Reshape array. Contour no down colums, point on crack along rows.
    end
elseif strcmpi(contourType,'none')
    %Take no action
else
    error('Unrecognised contourType. Should be either k, j, t, or none.')
end

%% Create output cell array
%If stepsToRead was not specified earlier, specify output for any step from
%0 to 1e6;
if isempty(stepsToRead)
    warning('stepsToRead not specified as input. Reading output for all available steps.')
    stepsToRead = [0:1e6];
end

%Build output cell array
if strcmpi(contourType,'k')
    outputLogical = ismember(kStep,stepsToRead);
    if exist('kContoursCell','var') && exist('k1Cell','var')
        outputArray = [kContoursCell(outputLogical);k1Cell(outputLogical);k2Cell(outputLogical);...
            k3Cell(outputLogical);kMTSDirCell(outputLogical);jFromKCell(outputLogical)];
    else
        warning('SIF output not read successfully from .dat file. Returning empty cell array.');
        outputArray = [{};{};{};{};{};{}];
    end
elseif strcmpi(contourType,'j')
    outputLogical = ismember(jStep,stepsToRead);
    if exist('jInt1Cell','var') && exist('jInt2Cell','var')
        outputArray = [jInt1Cell(outputLogical);jInt2Cell(outputLogical)];        %Returns cell array
    else
        warning('J-integral output not read successfully from .dat file. Returning empty cell array.');
        outputArray = [{};{}];
    end
elseif strcmpi(contourType,'t')
    outputLogical = ismember(tStep,stepsToRead);
    if exist('tStress1Cell','var') && exist('tStress2Cell','var')
        outputArray = [tStress1Cell(outputLogical);tStress2Cell(outputLogical)];          %Returns cell array
    else
        warning('T-stress output not read successfully from .dat file. Returning empty cell array.');
        outputArray = [{};{}];
    end
elseif strcmpi(contourType,'none')
    outputArray = [{};{}];  %This enables you to ask for no contour integral output, eg. if you want to get stepIncTimes alone.
end

%Optional output argument: Step, increment and time information.
if nargout == 2
    if exist('stepIncTimes','var')
        varargout{1} = stepIncTimes;
    else
        warning('Unable to extract step, increment and time information from the .dat file. Cannot return stepIncTimes array.')
    end
elseif nargout > 2
    error('Too many output arguments requested.')
end

end

