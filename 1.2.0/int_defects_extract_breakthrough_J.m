function [ ] = int_defects_extract_breakthrough_J( runName, modelNos, subsurfaceFlag, incNosOrTol, NFSDNodes, varargin )
%int_defects_extract_breakthrough_J.m
%Harry Coules 2018
%
%DESCRIPTION
%This function extracts the J-integral at various points on the crack tip
%line from models of embedded elliptical and surface semi-elliptical cracks
%in a plate of elastic-plastic material. The J-integral results are taken
%at the increment at which plastic breakthrough from the buried crack to
%the surface occurs. This is defined by a 'plastic tolerance', i.e. a
%von Mises equavalent plastic strain strain much be exceeded througout a
%path from the crack to the surface.
%
%INPUT ARGUMENTS
% runName - Name of the the int_defects run. There must be a folder called
%   abaqus_param_files_*runName* containing subdirectories which have
%   individual models in them.
% modelNos - A vector of model numbers from which J-integral results should
%   be extracted.
% subsurfaceFlag - A logical indicating whether the models are for
%   buried cracks (true) or surface semi-elliptical cracks (false). If the
%   cracks are on the surface, increment numbers for equivalent buried
%   cracks must be given as the third input argument.
% incNosOrTol - Either a (vector) list of increment numbers to extract
%   results from, or a (scalar) plastic strain tolerance. If it is a
%   plastic strain tolerance, int_defects_plastic_breakthrough_parametric
%   will be run to determine the plastic breakthrough increments.
% NFSDNodes - A vector of two node numbers. If subsurfaceFlag = true then
%   this should be the node number of the Nearest and Furthest point on
%   the crack tip line relative to the y = 0 plane (in that order). If
%   subsurfaceFlag = false, this should be the node number of a crack tip
%   line/Surface intersection node and of the Deepest point on the crack
%   tip line (in that order). This is used when reporting J-integral
%   results for specific points on the crack tip line.
%
%   *OPTIONAL*
% optionalSaveNameStr - A string which is included in the name of the saved
%   results file to distinguish it from other similar .mat files.
%
%OUTPUT ARGUMENTS
%   *none*
% - An output file containing the J-integral results will be written
%   (int_defects_extract_breakthrough_J_results.mat).
%
%NOTES
%   1. This function calls int_defects_extract_breakthrough_J.m three times
%   for every model being examined.
%   2. See int_defects_extract_breakthrough_J for the format of J-integral
%   results reported in outputArrayNear/Far/Surf/Deep/All.
%   3. The tolerance used for identifying contour-dependence in the contour
%   integral results is hard-coded in this function to 1e-2. See
%   int_defects_read_J_at_inc_parametric.m for more information.
%   4. For cracks woth 80 elements along the length of each semi-ellipse, 
%
%% Optional input argument
if nargin == 5
    optionalSaveNameStr = '';
elseif nargin == 6
    optionalSaveNameStr = varargin{1};  %Optional string to distinguish the saved results file.
else
    error('Unexpected number of input arguments.')
end    

%% Process final input argument. This is either a vector list of increment numbers or a scalar plastic breakthrough tolerance.
if length(incNosOrTol)>1
    incNosKnownFlag = true;
    incNos = incNosOrTol;       %If the final input argument has a length greater than 1, is must be the increment numbers.
elseif length(incNosOrTol)==1
    incNosKnownFlag = false;
    tol = incNosOrTol;          %If it has a length of exactly 1, is must be the plastic breakthrough tolerance.
else
    error('Final input argument has unexpected length')
end

%% Determine incNos
folderName = strcat('abaqus_param_files_',runName);
cd(folderName);

%Determine the plastic breakthrough increment for each model (if this is
%not already known)
if ~incNosKnownFlag
    [incNos,~,~,~] = int_defects_plastic_breakthrough_parametric( 0, tol, 'IntCrackJob1', modelNos );
end

%% Get J-integral results
contourTol = 1e-2;  %Note that a hard-coded contour-dependence tolerance of 1% is used.
if subsurfaceFlag
    [ outputArrayNear ] = int_defects_read_J_at_inc_parametric( modelNos, incNos, NFSDNodes(1), contourTol);    %Typ. 121
    [ outputArrayFar ] = int_defects_read_J_at_inc_parametric( modelNos, incNos, NFSDNodes(2), contourTol);     %Typ. 41
    [ outputArrayMax ] = int_defects_read_J_at_inc_parametric( modelNos, incNos, NaN, contourTol);  %Max. J anywhere on the crack tip line
    outputArrayAll = [outputArrayNear(:,[2:4]),outputArrayFar(:,[2:4]),outputArrayMax(:,[2:4])];
else
    [ outputArraySurf ] = int_defects_read_J_at_inc_parametric( modelNos, incNos, NFSDNodes(1), contourTol);    %Typ. 1
    [ outputArrayDeep ] = int_defects_read_J_at_inc_parametric( modelNos, incNos, NFSDNodes(2), contourTol);    %Typ. 41
    [ outputArrayMax ] = int_defects_read_J_at_inc_parametric( modelNos, incNos, NaN, contourTol);
    outputArrayAll = [outputArraySurf(:,[2:4]),outputArrayDeep(:,[2:4]),outputArrayMax(:,[2:4])];
end

cd ..

if isempty(optionalSaveNameStr)
    fileName = strcat('J_results_',runName,'.mat');
else
    fileName = strcat('J_results_',runName,optionalSaveNameStr,'.mat');     %Optional string to distinguish the saved results file.
end
save(fileName);

end