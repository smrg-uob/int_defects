function [] = int_defects_plate2cyl_auto( filename, curlAxis, pipeInnerRadius, radiusShiftFlag, varargin)
%int_defects_plate2cyl_auto.m
%Harry Coules 2016
%
%DESCRIPTION
%This function is used to generate an Abaqus input (.inp) file for a pipe
%containing a crack. It uses an input file generated for a plate
%containing an equivalent crack, and then applies a transformation so that
%this is curled around into a hollow cylinder.
%
%INPUT ARGUMENTS
%   filename - Name of the file to be modified as a string.
%   curlAxis - String defining the model axis about which the mesh should
%       be curled. Can be 'x', 'y' or 'z'.
%   pipeInnerRadius - Internal radius of the cylinder. Defaults to 0.5225
%       if this input argument is not supplied.
%   radiusShiftFlag - Logical indicating whether the initial coordinate
%       system needs to be shifted to account for pipe inner radius. Should be
%       set to true if the initial coordinate system origin is located at the
%       plate's surface, i.e. when using models generated with int_defects.
%*OPTIONAL*
%   antiCurlFlag - Logical indicating that the pipe should be curled in the
%       anticlockwise direction (so that Crack A is on the outer diameter).
%       If this is not supplied, the default is false.
%   pipeWallThickness - Cylinder wall thickness. This must be supplied if
%       antiCurlFlag = true.
%
%OUTPUT ARGUMENTS
%   *none*
%
%NOTES
% - When multiple optional input arguments are given, they must be given in
% the order: antiCurlFlag, pipeWallThickness.
%
%% Optional input argument
if isempty(varargin)
    antiCurlFlag = false;
    pipeWallThickness = NaN;
elseif length(varargin) == 2
    antiCurlFlag = varargin{1};
    pipeWallThickness = varargin{2};
else
    error('Unexpected number of input arguments.')
end

%% Scan the existing .inp file to get the node definitions
tic
disp('TRANSFORMING PLATE GEOMETRY TO PIPE GEOMETRY')
disp('Scanning existing .inp file...');
inpScanInputStruct.filename = filename;
[inpStruct] = read_inp(inpScanInputStruct);

%Check that nodes are uniquely numbered
if length(inpStruct.nodes(:,1)) ~= length(unique(inpStruct.nodes(:,1)))
    error('Duplicate nodes have been detected in the .inp file. Cannot proceed.')
end

%% Calculate new coordinates
disp('Calculating new nodal coordinates...');
if strcmpi(curlAxis,'x')
    rCol = 3;
    phiCol = 4;
elseif strcmpi(curlAxis,'y')
    rCol = 4;
    phiCol = 2;
elseif strcmpi(curlAxis,'z')
    rCol = 3;
    phiCol = 2;
end

%Calculate new nodal coords
for k1 = 1:size(inpStruct.nodes,1)   
    if antiCurlFlag
        if radiusShiftFlag
            r = inpStruct.nodes(k1,rCol) - pipeInnerRadius - pipeWallThickness;
            phi = -inpStruct.nodes(k1,phiCol)./pipeInnerRadius;
        else
            r = inpStruct.nodes(k1,rCol);
            phi = -inpStruct.nodes(k1,phiCol)./pipeInnerRadius;
        end
    else
        if radiusShiftFlag
            r = inpStruct.nodes(k1,rCol) + pipeInnerRadius;
            phi = inpStruct.nodes(k1,phiCol)./pipeInnerRadius;
        else
            r = inpStruct.nodes(k1,rCol);
            phi = inpStruct.nodes(k1,phiCol)./pipeInnerRadius;
        end
    end
    inpStruct.nodes(k1,rCol) = r.*cos(phi);
    inpStruct.nodes(k1,phiCol) = r.*sin(phi);
end

%% Determine the start and end of the node definitions in the original .inp file
% Read in the lines of the .inp file
fid = fopen(filename,'r','n','UTF-8');
C = textscan(fid, '%s', Inf, 'Delimiter','\n'); C = C{1};
fclose(fid);

%Determine the start and end lines of the node definitions
nodeDefEndIndexFlag = false;
k1 = 1;
while ~nodeDefEndIndexFlag && k1 <= length(C)
    if strncmpi(char(C(k1)),'*Node',5);
        nodeDefStartIndex = k1;
    elseif strncmpi(char(C(k1)),'*Element',8);
        nodeDefEndIndex = k1-1;
        nodeDefEndIndexFlag = true;
    elseif k1 == length(C)
        error('Reached the end of the input file without finding the end of the node definitions.')
    end
    k1 = k1 + 1;
end

%% Write new .inp file
%Rename the old .inp file
disp('Renaming old (plate geometry) .inp file...');
[token,remain] = strtok(filename,'.');
movefile(filename,[token,'_plate',remain]);

%Write new
disp('Writing new (cylinder geometry) .inp file...');
inpWriteStruct.NodeDefs.stringLines = {'*Node\n'};
inpWriteStruct.NodeDefs.formatSpec = '%d, %E, %E, %E\n';
inpWriteStruct.NodeDefs.data = inpStruct.nodes;

fid = fopen(filename,'w');
for k1 = 1:length(C)
    if k1 < nodeDefStartIndex
        fprintf(fid,'%s\n',char(C(k1)));
    elseif k1 == nodeDefStartIndex
        fprintf(fid,inpWriteStruct.NodeDefs.stringLines{1});
    elseif k1 == nodeDefStartIndex + 1
        fprintf(fid,inpWriteStruct.NodeDefs.formatSpec,inpWriteStruct.NodeDefs.data');
    elseif k1 > nodeDefEndIndex
        fprintf(fid,'%s\n',char(C(k1)));
    end
end
fclose(fid);

disp(['int_defects_plate2cyl_auto terminated successfully after ',num2str(toc,'%4.2f'),' seconds']);

end