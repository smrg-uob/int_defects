function [ incNo, incTime, modelMaxInc, modelMaxTime ] = int_defects_plastic_breakthrough( startPoints, endPlaneY, ePl, jobName )
%int_defects_plastic_breakthrough.m
%Harry Coules 2018
%
%DESCRIPTION
%This function determines the model increment at which "plastic breaktrough"
%occurs in an Abaqus finite element model. Plastic breakthrough is defined
%as the appearance of a unbroken path of elements, each with a plastic
%strain value greater than a specified tolerance, which connects a specified
%location (or any one of a set of locations) with a specified plane. This
%is useful in fracture mechanics analyses where it is necessary to
%determine the time at which the strain across an un-cracked ligament in a
%structure reaches an unacceptable value.
%
%INPUT ARGUMENTS
%   startPoints - n-by-3 array specifing n start locations (usually nodes
%       on the crack tip line) in 3D space.
%   endPlaneY - y-coordinate of the plane to which plastic break-through is
%       being determined. It is assumed that this lies parallel to the x-z
%       plane. For a pipe model, this should be the y-coordinate of the
%       intersection of the pipe outer face with the crack plane.
%   ePl - The plastic strain tolerance. The plastic strain must be greater
%       than this value throughout a path from the starting point to the
%       ending surface.
%   jobName - Name of the Abaqus job.
%
%OUTPUT ARGUMENTS
%   incNo - The increment number in which plastic breakthrough first
%       occurs, given the tolerance specified by ePl.
%   incTime - The (step) time at which plastic breakthrough first occurs.
%   modelMaxInc - The last increment in the model.
%   modelMaxTime - The (step) time at which the last increment occurs.
%
%USAGE
% 1. Navigate to the directory containing the Abaqus .dat and .inp files.
% 2. Run this function.
%
%NOTES
% 1. This function assumes certain things about the Abaqus model. For
%   example that there are two types of element and that these are linear
%   brick (8-noded) and quadratic tet (10-noded). It also assumes that the
%   ending surface is made of tets only. It is assumed that the loading is
%   monotonic, so that if plastic breakthrough is not present in the final
%   increment of the model, then it is not present in any increment. This
%   function is designed to be compatible with models automatically
%   generated using the int_defects package, but its applicability to other
%   models is untested.
% 2. If plastic breakthrough is determined not to have occurred, this
%   function returns incNo = NaN.
% 3. When multiple startPoints are given, a plastic path from any of them
%   to the specified plane is allowed. Therefore the incNo returned by this
%   function indicates the first increment in which any one of them is
%   connected to the plane.
%
%% Read the input file
tic;
disp('Reading and extracting information from input file. This may take some time...');

%If the geometry has been transformed from plate to cylindrical geometry,
%read the nodal coordinates from the plate .inp file.
if ~isempty(dir('*_plate.inp'))
    inpScanInputStruct.filename = strcat(jobName,'_plate.inp');
else
    inpScanInputStruct.filename = strcat(jobName,'.inp');
end
[inpStruct] = read_inp(inpScanInputStruct);         %Scan input file

%If there is a Reference Point node (for load application) remove it from the table
node1Idx = inpStruct.nodes(:,1) == 1;               %All nodes with node no. = 1
if sum(node1Idx)>1
    [~,idx] = max(node1Idx.*inpStruct.nodes(:,4));  %The Node 1 with the maximum z coordinate
    inpStruct.nodes(idx,:) = [];                    %Remove it.
end

%Determine starting node no.s
noStartPoints = size(startPoints,1);                                %Number of start points which plastic connectivity will be tested for.
[uniqueNodes,idx] = unique(inpStruct.nodes(:,2:end),'rows');        %Uniquely-located nodes (Note: there may be duplicate notes at tie constraints)
delaunayTriNodes = delaunayTriangulation(uniqueNodes);              %Delaunay triangulation
startNodes = nearestNeighbor(delaunayTriNodes,startPoints(:,1),startPoints(:,2),startPoints(:,3));    %Use nearest-neighbour nodes to starting locations as start nodes
startNodes = idx(startNodes);                                                                         %Index in set of ALL nodes from index in set of uniquely-located nodes

%Determine ending node no.s
idx = inpStruct.nodes(:,3) == endPlaneY;
endNodes = inpStruct.nodes(idx,1);

%Create a table of all elements with node assignments. Pad with NaNs.
nanFill = NaN*ones(size(inpStruct.elems{1,1},1),2);
elems = [[inpStruct.elems{1,1},nanFill];inpStruct.elems{1,2}];  %Note this is hard-coded to expect linear brick then quadratic tet elements.
elems = sortrows(elems);
if elems(end,1) ~= size(elems,1)
    error('Unexpected number of rows in elements table.');
end

%Elements that lie on the ending plane
endElems = [];
for k1 = 1:size(elems,1)
    C = intersect(endNodes,elems(k1,2:end));
    if length(C)==6
        endElems = [endElems;elems(k1,1)];
    elseif length(C)>6
        error('An element has an unexpected number of nodes on the ending plane.')
    end
end

disp('... finished extracting information from input file.');

%% Check what the maximum increment that the model actually completed was
datScanInputStruct.filename = strcat(jobName,'.dat');
[ ~, stepIncTimes ] = read_dat_contour_integral( datScanInputStruct,'none',1);
modelMaxInc = stepIncTimes(end,2);
modelMaxTime = stepIncTimes(end,6); %Note max. step time is used

%% Determine the plastic breakthrough increment
%Loop over increments, checking the plastic breakthrough condition in each.
critSatisfied = NaN*ones(1,modelMaxInc);    %Vector containing plastic breakthrough criterion result for each increment
incNoCheck = floor(modelMaxInc/2);          %First increment to try
incStep = floor(modelMaxInc/2);             %Initial step size
incFoundFlag = false;                       %Stopping criterion
while incFoundFlag == false
    disp(['Checking for plastic breakthrough in increment no. ',num2str(incNoCheck),'...']);
    
    %Scan the dat file
    datScanInputStruct.filename = strcat(jobName,'.dat');
    datScanInputStruct.stepsToRead = 1;
    datScanInputStruct.incsToRead = incNoCheck;
    [datStruct] = read_dat(datScanInputStruct);
    
    %Determine whether the plastic breakthrough criterion has been
    %satisfied on this increment.
    breakthroughFlag = false;
    if ~isempty(datStruct.peeq)
        peeq = datStruct.peeq(datStruct.peeq(:,2)>ePl,:);           %Elements containing peeq above the required amount
        nodeConnect = elems(ismember(elems(:,1),peeq(:,1)),2:end);  %Nodes in these elements
        
        if ~isempty(peeq)
            %Form an unweighted graph of plastically-connected nodes
            nodeConnect2 = zeros(size(nodeConnect,1)*9,2);  %Preallocate nodeConnect2
            for k1 = 1:size(nodeConnect,1)
                nodeConnectBlock = [nodeConnect(k1,1)*ones(9,1),nodeConnect(k1,2:10)']; %For each element, the first node is connected to all others
                nodeConnect2([(k1*9)-8:(k1*9)],:) = nodeConnectBlock;
            end
            nodeConnect = nodeConnect2;
            clear nodeConnect2
            
            nodeConnect(isnan(nodeConnect(:,1)),:) = [];  %Reshape and remove NaNs
            nodeConnect(isnan(nodeConnect(:,2)),:) = [];
            nodeConnect = unique(sortrows(nodeConnect),'rows'); %Re-order and find unique rows
            
            duplicateFlag = ismember(nodeConnect,fliplr(nodeConnect),'rows');   %Remove duplicate graph edges
            nodeConnect(duplicateFlag,:) = [];
            
            G = graph(nodeConnect(:,1),nodeConnect(:,2));   %Create unweighted graph
            
            %Now test the plastic breakthrough criterion - for each start point.
            for k1 = 1:noStartPoints
                %Solve the graph connectivity
                connNodes = dfsearch(G,startNodes(k1));                   %All nodes plastically connected to the start point (depth-first search)
                connSurfNodes = connNodes(ismember(connNodes,endNodes));  %Connected surface nodes
                
                %If we can connect to 6 surface nodes of the same element, the
                %criterion is satisfied. Note this assumes that the surface
                %element is a quadratic tet.
                for k2 = 1:length(endElems)
                    endElemNodes = elems(endElems(k2),2:end);
                    C = intersect(connSurfNodes,endElemNodes);
                    if length(C)==6
                        breakthroughFlag = true;
                    end
                end
            end
            
            %Add info on whether or not breakthrough was detected in this
            %increment to the output vector.
            if breakthroughFlag
                critSatisfied(incNoCheck) = 1;                  %Set flag for criterion satisfied in this increment...
                incStep = -min([ceil(abs(incStep)/2),incNoCheck]);   %... and set incStep so that we'll check an earlier increment next iteration.
            else
                critSatisfied(incNoCheck) = 0;                              %Set flag for criterion not satisfied in this increment...
                incStep = min([ceil(abs(incStep)/2),modelMaxInc-incNoCheck]);    %... and set incStep so that we'll check an later increment next iteration.
            end
        else
            critSatisfied(incNoCheck) = 0;
            incStep = min([ceil(abs(incStep)/2),modelMaxInc-incNoCheck]);
        end
    else
        critSatisfied(incNoCheck) = 0;
        incStep = min([ceil(abs(incStep)/2),modelMaxInc-incNoCheck]);
    end
    
    %Check if a turning point between non-satisfied and satisfied has
    %appeared in critSatisfied
    if ~isempty(strfind(critSatisfied,[0,1]))
        incFoundFlag = true;
        incNo = strfind(critSatisfied,[0,1])+1;
    elseif incNoCheck == modelMaxInc
        incFoundFlag = true;
        incNo = NaN;    %If we reach the last increment in the model without the breakthrough criterion having been satisfied, return incNo = NaN.
    end
    
    %Calculate the next increment to check
    incNoCheck = incNoCheck + incStep;
end

%Determine the step time at which the breakthrough increment occurs.
if isnan(incNo)
    incTime = NaN;
    disp(['Plastic breakthrough does not occur.'])
else
    incTime = stepIncTimes(incNo,5);
    disp(['Breakthrough occurs at increment no. ',num2str(incNo),' (increment time: ',num2str(incTime),').']);
end
disp(['Time taken to determine breakthrough increment for this model: ',num2str(toc,'%4.2f'),' seconds']);
disp(' ');

end