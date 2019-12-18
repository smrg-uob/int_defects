function [ incNos, incTimes, modelMaxIncs, modelMaxTimes ] = int_defects_plastic_breakthrough_parametric( endPlaneY, ePl, jobName, modelNos )
%int_defects_plastic_breakthrough_parametric.m
%Harry Coules 2018
%
%DESCRIPTION
%This function loops through a set of elastic-plastic models of a single
%embedded crack, determining the increment at which plastic breakthrough
%occurs to a given surface.
%
%INPUT ARGUMENTS
%   endPlaneY - y-coordinate of the plane to which plastic break-through is
%       being determined. It is assumed that this lies parallel to the x-z
%       plane, and is the same for all models.
%   ePl - The plastic strain tolerance. The plastic strain must be greater
%       than this value throughout a path from the starting point to the
%       ending surface. An arbitrarily-low value (eg. 1e-6) is recommended
%       in most cases.
%   jobName - Name of the Abaqus job.
%   modelNos - A list of model numbers to determine the breakthrough
%       increment for.
%
%OUTPUT ARGUMENTS
%   incNos - a vector containing the model increments at which plastic
%       breakthrough is judged to have occurred.
%   incTimes - a vector containing the (step) times at which plastic
%       breakthrough is judged to have occurred.
%   modelMaxIncs - a vector containing the last increments of each model.
%   modelMaxTimes - a vector containing the (step) times at which the last
%       increments in the model occur.
%
%USAGE
% 1. Navigate to the abaqus_param_files directory containing the model
%   sub-directories.
% 2. Run this function.
%
%NOTES
% 1. This function defines the startpoints for the plastic breakthrough
%   check using the extreme points in y of the cracks. It determines what
%   crack are present, and their geoemtry, from
%   paramLogStruct(modelNo).modelParams.
% 2. As a result of the many assumptions that this function and
%   int_defects_plastic_breakthrough make about the model and model
%   geometry, it is considered unlikely that it will work for any FE model
%   output other than that created using the int_defects software package.
% 3. When multiple cracks are present, this function determines the first
%   model increment when plastic breakthrough occurs from ANY of the cracks
%   to the back face.
%
%%
tic;

%Load the workspace_dump file so that we can extract the starting point
%locations.
load('workspace_dump.mat','paramLogStruct');

%Loop through model sub-directories
incNos = [];
incTimes = [];
modelMaxIncs = [];
modelMaxTimes = [];
for k1 = 1:length(modelNos)
    modelNo = modelNos(k1);
    disp(' ');
    disp(['Determining plastic breakthrough increment for model #',num2str(modelNo)]);
    disp(['This is ',num2str(k1),' out of ',num2str(length(modelNos)),' in the sequence to scan.']);
    
    %Initially set breakthrough increment etc. as NaN. This will be output if the model did not run.
    %Note: jobExitStatus == 2 indicates that it was not necessary to run the model.
    %Note: We would expect jobExitStatus == 0 for most models. Generally, limit load models will fail when the deformation becomes unbounded.
    incNo = NaN;
    incTime = NaN;
    modelMaxInc = NaN;
    modelMaxTime = NaN;

    if ~(paramLogStruct(modelNo).jobExitStatus==2)    %i.e. If a model with this number has run.
        %Calculate the location(s) of the starting point(s) on the crack tip line.
        startPoint = [];
        if paramLogStruct(modelNo).modelParams.singleCrackFlag
            %Single crack
            if paramLogStruct(modelNo).modelParams.crackPositions.subsurfaceAFlag   %Single embedded crack
                %Furthest point from y=0 plane
                startPoint(1,1) = paramLogStruct(modelNo).modelParams.crackPositions.xA;
                startPoint(1,2) = paramLogStruct(modelNo).modelParams.crackPositions.yA + paramLogStruct(modelNo).modelParams.crackSizes.aA1;
                startPoint(1,3) = 0;
                %Nearest point to y=0 plane
                startPoint(2,1) = paramLogStruct(modelNo).modelParams.crackPositions.xA;
                startPoint(2,2) = paramLogStruct(modelNo).modelParams.crackPositions.yA - paramLogStruct(modelNo).modelParams.crackSizes.aA1;
                startPoint(2,3) = 0;
            else   %Single surface crack
                %Furthest point from y=0 plane
                startPoint(1,1) = paramLogStruct(modelNo).modelParams.crackPositions.xA;
                startPoint(1,2) = paramLogStruct(modelNo).modelParams.crackPositions.yA + paramLogStruct(modelNo).modelParams.crackSizes.aA1;
                startPoint(1,3) = 0;
            end
        else
            %Two cracks
            if paramLogStruct(modelNo).modelParams.crackPositions.oppositeFlag
                if paramLogStruct(modelNo).modelParams.crackPositions.subsurfaceBFlag   %Two cracks - one surface, one embedded.
                    %Furthest point from y=0 plane - Crack A
                    startPoint(1,1) = paramLogStruct(modelNo).modelParams.crackPositions.xA;
                    startPoint(1,2) = paramLogStruct(modelNo).modelParams.crackPositions.yA + paramLogStruct(modelNo).modelParams.crackSizes.aA1;
                    startPoint(1,3) = 0;
                    %Furthest point from y=0 plane - Crack B
                    startPoint(2,1) = paramLogStruct(modelNo).modelParams.crackPositions.xB;
                    startPoint(2,2) = paramLogStruct(modelNo).modelParams.crackPositions.yB + paramLogStruct(modelNo).modelParams.crackSizes.aB1;
                    startPoint(2,3) = 0;
                    %Nearest point to y=0 plane - Crack B
                    startPoint(3,1) = paramLogStruct(modelNo).modelParams.crackPositions.xB;
                    startPoint(3,2) = paramLogStruct(modelNo).modelParams.crackPositions.yB - paramLogStruct(modelNo).modelParams.crackSizes.aB1;
                    startPoint(3,3) = 0;
                else
                    error('paramLogStruct indicates that there are two surface cracks on opposite surfaces. This function cannot be used for this geometry.');
                end
            else   %Two surface cracks
                %Furthest point from y=0 plane - Crack A
                startPoint(1,1) = paramLogStruct(modelNo).modelParams.crackPositions.xA;
                startPoint(1,2) = paramLogStruct(modelNo).modelParams.crackPositions.yA + paramLogStruct(modelNo).modelParams.crackSizes.aA1;
                startPoint(1,3) = 0;
                %Furthest point from y=0 plane - Crack B
                startPoint(2,1) = paramLogStruct(modelNo).modelParams.crackPositions.xB;
                startPoint(2,2) = paramLogStruct(modelNo).modelParams.crackPositions.yB + paramLogStruct(modelNo).modelParams.crackSizes.aB1;
                startPoint(2,3) = 0;
            end            
        end
        
        cd(num2str(modelNo,'%06i')); %Navigate to model directory
        [incNo, incTime, modelMaxInc, modelMaxTime] = int_defects_plastic_breakthrough( startPoint, endPlaneY, ePl, jobName );    %Determine breakthrough increment.
        cd ..   %Return to main directory
    end
    
    incNos = [incNos;incNo];
    incTimes = [incTimes;incTime];
    modelMaxIncs = [modelMaxIncs;modelMaxInc];
    modelMaxTimes = [modelMaxTimes;modelMaxTime];
end

end