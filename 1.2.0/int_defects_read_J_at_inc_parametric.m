function [ outputArray ] = int_defects_read_J_at_inc_parametric( modelNo, maxInc, nodeNo, varargin)
%int_defects_read_J_at_inc_parametric.m
%Harry Coules 2018
%
%DESCRIPTION
%This function is used for extracting J-integral results from a series of
%results from Abaqus models of cracked bodies. It loops through models and
%extracts contour integral output for a single specified location, which
%must be the same in each model.
%
%INPUT ARGUMENTS
% modelNo - Ordered vector of model numbers.
% maxInc - List of increments at which the criterion for extracting the
%   J-integral has been satisfied. Should be a vector of the same length as
%   modelNo.
% nodeNo - The node number for which J-integral results are required. For
%   an embedded crack with 161 nodes on the crack tip line, this is
%   typically 41 (Jfar) or 121 (Jnear). nodeNo can be given as NaN, in
%   which case the maximum J at any point on the crack tip line will be
%   returned.
%
% *OPTIONAL*
% tol - Relative tolerance on the J-integral path-dependence. This
%   specifies the maximum relative difference between contour integral
%   results from the last and penultimate contours. Default is 0.01.
%
%OUTPUT ARGUMENTS
% outputArray - An n-by-5 array containing (in columns):
%   1. (Best-estimate) remotely-applied stress.
%   2. (Best-estimate) J-integral at requested location.
%   3. Lower limit of J-integral at requested location (from previous inc).
%   4. Upper limit of J-integral at requested location (from current inc).
%   5. A logical which flags J-integral path-dependence.
%
%USAGE
% 1. Navigate into the abaqus_param_files directory containg the model
% results which will be interrogated.
% 2. Run int_defects_read_J_at_inc_parametric.
%
%NOTES
%   1. This function is called by int_defects_extract_breakthrough_J.m.
%
%% Optional input argument
if isempty(varargin)
    tol = 0.01; %The default tolerance is 0.01.
elseif length(varargin) == 1
    tol = varargin{1};
else
    error('Too many input arguments.')
end

%% Loop through models, determining the relevant J result for each model.
outputArray = zeros(length(modelNo),5); %Preallocate output array

%Loop over models
for k1 = 1:length(modelNo)
    if isnan(maxInc(k1))                        %In the case where a model run has not been attempted, or the model has not reached plastic breakthorugh...
         outputArray(k1,:) = NaN*ones(1,5);     %... add a line of NaNs to the output array.
    else
        cd(num2str(modelNo(k1),'%06i'));                                                                             %Otherwise, navigate to model directory
        
        [ J, JLims, remoteStress, pathDepWarnFlag ] = int_defects_read_J_at_inc( maxInc(k1), nodeNo, tol);  %Read J-integral information
        outputArray(k1,:) = [remoteStress, J, JLims, pathDepWarnFlag ];                                     %Add result to output array
        
        cd ..   %Return to main directory
    end
end

end