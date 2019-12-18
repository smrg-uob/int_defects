function [ outputStruct ] = read_inp( inputStruct )
%read_inp.m
%Harry Coules 2013
%   Updates:
%       - 04/2019: Improved array preallocation for better performance.
%
%DESCRIPTION
%This function scans an Abaqus/Standard .inp file and extracts information
%about the nodes and elements of the solid model, i.e. their locations and
%connectivity.
%
%INPUT ARGUMENTS (Input structure may contain:)
%   inputStruct.filename - Name of the .inp file to be read.
%   inputStruct.jobName - Name of Abaqus job.
%       NOTE: It is necessary to specify one of the above fields.
%   inputStruct.preallocateLength - Length for which node and element
%       arrays will be preallocated. If the number of nodes in the input
%       file exceeded preallocateLength, severe performance loss would
%       occur, so this is prevented in the code (it will throw an error).
%       Therefore preallocateLength should be set to a larger number
%       than the number of nodes in the model which will be read. The
%       default (if this field is not supplied) is 1e6.
%
%OUTPUT ARGUMENTS (Output structure can contain:)
%   outputStruct.nodes - Node no.s and nodal coordinates
%   outputStruct.elems - Element no.s and corresponding node no.s for each
%       element.
%
%NOTES:
%1. The nodal coordinates are given in the global coordinate system of the
%   Abaqus assembly. Depending on how you've set up the model, this may or
%   may not coincide with the coordinate system of the part.
%2. Currently this function only works for .inp files in which the assemply
%   contains only a single part.
%3. This function only works for parts in which a single element type is
%   used.
%4. Care should be taken to ensure that inputStruct.preallocateLength is
%   set to an appropriate value, if needed. Note that this is limited by
%   the available system memory.
%
%% Check input structure
%If filename and jobName are empty in input structure, give an error.
if isfield(inputStruct,'filename')
    filename = inputStruct.filename;
elseif isfield(inputStruct,'jobName')
    filename = strcat(inputStruct.jobName,'.inp');
else
    error('No filename or job name defined in input structure to specify the .inp file.');
end

%Array preallocation length
if isfield(inputStruct,'preallocateLength')
    preallocateLength = inputStruct.preallocateLength;
else
    preallocateLength = 1e6;
end

%% Read in the lines of the .inp file
fid = fopen(filename,'r','n','UTF-8');
C = textscan(fid, '%s', Inf, 'Delimiter','\n'); C = C{1};
fclose(fid);

%% Scan the .inp file
nodes = NaN*ones(preallocateLength,4);
elems = NaN*ones(preallocateLength,32);
inputType = '';
elemDefTableNo = 0;
elemDefFirstLineFlag = true;
nodeCount = 0;
elemCount = 0;
for k = 1:length(C)
    if nodeCount >= preallocateLength || elemCount >= preallocateLength
        error(['No. of nodes or elements in the model exceeds preallocated array length of ',num2str(preallocateLength)]);
    end
    
    %Determine if you're at the start or end of the nodes/elements
    %definition section
    if strncmpi(char(C(k)),'*Node',5);
        inputType = 'NODE';
        nodeCount = 0;
    
    elseif strncmpi(char(C(k)),'*Element',8) && ~strncmpi(char(C(k)),'*Element Output',15); 
        inputType = 'ELEMENT';
        elemDefTableNo = elemDefTableNo + 1;
        elemCount = 0;
        if elemDefTableNo == 2
            %If a second table of element set definitions is detected,
            %we need to use cell array rather than table output for elems
            elemsTemp = elems;
            clear elems
            elems{1} = elemsTemp;
            clear elemsTemp
            elems{elemDefTableNo} = NaN*ones(preallocateLength,32);
        elseif elemDefTableNo > 2
            elems{elemDefTableNo} = NaN*ones(preallocateLength,32);
        end
        
        %Now check what type of element is being used, so that we can tell
        %how many nodes it will take to define each one.
        if isscalar(strfind(char(C(k)), 'CAX4R'));  %Axisymmetric, 4-node bilinear, reduced integration with hourglass control
            elemStrLength = 5;
        elseif isscalar(strfind(char(C(k)), 'CPS4R'));  %Plane stress, 4-node bilinear, reduced integration with hourglass control
            elemStrLength = 5;
        elseif isscalar(strfind(char(C(k)), 'CPE4R'));  %Plane strain, 4-node bilinear, reduced integration with hourglass control
            elemStrLength = 5;
        elseif isscalar(strfind(char(C(k)), 'CPE8R'));  %Plane strain, 8-node biquadratic
            elemStrLength = 9;
        elseif isscalar(strfind(char(C(k)), 'CPS8R'));  %Plane stress, 8-node biquadratic
            elemStrLength = 9;
        elseif isscalar(strfind(char(C(k)), 'C3D8'));  %3D, 8-node linear brick, reduced integration with hourglass control
            elemStrLength = 9;
        elseif isscalar(strfind(char(C(k)), 'C3D4'));   %3D, 4-node linear
            elemStrLength = 5;
        elseif isscalar(strfind(char(C(k)), 'C3D6'));   %3D, 6-node linear triangular prism
            elemStrLength = 7;
        elseif isscalar(strfind(char(C(k)), 'C3D10'));  %3D, 10-node quadratic tetrahedron
            elemStrLength = 11;
        elseif isscalar(strfind(char(C(k)), 'C3D15'));  %3D, 15-node quadratic triangular prism
            elemStrLength = 16;
        elseif isscalar(strfind(char(C(k)), 'C3D20'));  %3D, 20-node quadratic brick
            elemStrLength = 21;
        elseif isscalar(strfind(lower(char(C(k))), 'output'));  %In this case, you've reached the *Element Output line
            inputType = '';
        else
            error('The type of element used in the model is unfamiliar. Cannot extract element definitions from the .inp file.')
        end
        
    elseif strncmpi(char(C(k)),'*Instance',9);
        inputType = 'INSTANCE';
        
    elseif strncmpi(char(C(k)),'*Nset',5)||strncmpi(char(C(k)),'*Elset',6)||strncmpi(char(C(k)),'*End Assembly',13)||strncmpi(char(C(k)),'*End Instance',13);
        inputType = '';
        
    %Add node/element/instance definition to the matrix, if appropriate
    elseif strcmp(inputType,'NODE')    %Format for node definitions
        line = sscanf(char(C(k)),'%u, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f')';
        if length(line) == 3
            nodeCount = nodeCount+1;
            nodes(nodeCount,1:3) = line;   %Construct matrix of nodes & nodal coordinates
        elseif length(line) == 4
            nodeCount = nodeCount+1;
            nodes(nodeCount,1:4) = line;   %Construct matrix of nodes & nodal coordinates
        end
        
    elseif strcmp(inputType,'ELEMENT')    %Format for element definitions
        line = sscanf(char(C(k)),'%u, %u, %u, %u, %u, %u, %u, %u, %u, %u, %u, %u, %u, %u, %u, %u')';
        if elemStrLength <= 16
            if length(line) == elemStrLength;
                elemCount = elemCount+1;
                if elemDefTableNo == 1
                    elems(elemCount,1:elemStrLength) = line;   %Construct matrix of elements with corresponding nodes
                else
                    elems{elemDefTableNo}(elemCount,1:elemStrLength) = line;
                end
            end
        elseif elemStrLength <= 32
            if elemDefFirstLineFlag
                %First line of the element definition - add to data,
                %padding with NaNs.
                if length(line) == 16;
                    elemCount = elemCount+1;
                    if elemDefTableNo == 1
                        elems(elemCount,1:elemStrLength) = [line,NaN*zeros(1,elemStrLength-16)];   %Construct matrix of elements with corresponding nodes
                    else
                        elems{elemDefTableNo}(elemCount,1:elemStrLength) = [line,NaN*zeros(1,elemStrLength-16)];
                    end
                else
                    error('When an element definition definition runs to two lines, the first line should contain exactly 16 integers.')
                end
                elemDefFirstLineFlag = false;
            else
                %Second line of the element definition - add to data
                %structure where the NaNs are.
                if length(line) == elemStrLength-16;
                    if elemDefTableNo == 1
                        elems(elemCount,17:elemStrLength) = line;   %Construct matrix of elements with corresponding nodes
                    else
                        elems{elemDefTableNo}(elemCount,17:elemStrLength) = line;
                    end
                else
                    error('When an element definition definition runs to two lines, the second line should contain exactly elemStrLength-16 integers.')
                end
                elemDefFirstLineFlag = true;
            end
        else
            error('This function cannot read node/element data for elements with more than 31 nodes.')
        end
        
    elseif strcmp(inputType,'INSTANCE')    %Format for part instance position
        line = sscanf(char(C(k)),'%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f')';
        if length(line) == 3
            if ~exist('instancePos','var')  %First line defines a translation for the part instance
                instancePos = line;
            elseif ~exist('instanceRot','var')  %Second line (optional) defines a rotation
                instanceRot = line;
            else
                warning('Too many data lines for part instance positioning. Expected: 1 or 2 lines.')
                warning('The node & element coordinates determined by read_inp may be incorrect!')
            end
        end
    end
end

%% Remove any excess preallocated rows and columns from nodes and elems
nodes(:,~any(~isnan(nodes))) = [];
nodes(~any(~isnan(nodes),2),:) = [];

if isnumeric(elems)     %If elems is an array
    elems(:,~any(~isnan(elems))) = [];
    elems(~any(~isnan(elems),2),:) = [];
else
    for k1 = 1:length(elems)    %If elems is a cell array, loop over cells
        elems{k1}(:,~any(~isnan(elems{k1}))) = [];
        elems{k1}(~any(~isnan(elems{k1}),2),:) = [];
    end
end

%% Organise output into a structure
%Ensure all output data is sorted by element number
if elemDefTableNo >1
    for k1 = 1:length(elems)
        if ~isempty(elems{k1})
            elems{k1} = sortrows(elems{k1},1);
        end
    end
else
    if ~isempty(elems)
        elems = sortrows(elems,1);
    end
end

if ~isempty(nodes)
    nodes = sortrows(nodes,1);
end

%If the part instance is offset from the model centre, apply this
%translation to the nodal coordinates.
if exist('instancePos','var')
    nodes(:,2:4) = nodes(:,2:4)+ones(size(nodes,1),1)*instancePos;
end
%If there is a rotation, don't apply it (because this isn't necessary at
%the moment and I can't be arsed with the algebra).
if exist('instanceRot','var')
    warning('A part instance rotation has been detected in the .inp file, but this has not been applied to the returned nodal coordinates.');
end

%Construct a structure containing the extracted data
outputStruct.filename = filename;
outputStruct.nodes = nodes;
outputStruct.elems = elems;
outputStruct.elemDefTableNo = elemDefTableNo;

end