function [ J, JLims, modelTime, pathDepWarnFlag ] = int_defects_read_J_at_inc( maxInc, nodeNos, varargin)
%int_defects_read_J_at_inc.m
%Harry Coules 2018
%
%DESCRIPTION
%This function is used for determining J-integral values for at a specific
%location and point in time from a set of crack driving force models run
%using the int_defects package.
%
%INPUT ARGUMENTS
% maxInc - Model increment at which J-integral results will be extracted.
%   Should be a scalar integer.
% nodeNos - A vector of node numbers for which J-integral results are
%   required. If nodeNo is given as NaN, the function will determine the
%   location of the greatest J anywhere on the crack tip line and return
%   the result for that location.
% *OPTIONAL*
% tol - Relative tolerance on the J-integral path-dependence. This
%   specifies the maximum relative difference between contour integral
%   results from the last and penultimate contours. Default is 0.01.
%
%OUTPUT ARGUMENTS
% J - Best-estimate for the J-integral at the time at which the extraction
%   criterion was satisfied (i.e. at the mean of the J results at maxInc
%   and maxInc-1).
% JLims - Min. and max. estimates for the J-integral at the time when the
%   extraction criterion was specified. These bounds are simply the
%   J-integral results for the increments which bracket when the criterion
%   was satisfied.
% modelTime - Model time at the moment of the J-integral best-guess.
% pathDepWarnFlag - Logical to indicate that there may be path-dependece
%   issues with the contour inegral results.
%
%USAGE
% 1. Determine the model increment maxInc. This migh be done eg. by visual
%   inspection of an Abaqus .odb file to determine at which increment the
%   extraction criterion has been satisfied. Or through automatic detection
%   by int_defects_plastic_breakthrough.
% 2. In Matlab, navigate the the subdirectory which contains results for
%   the model in question.
% 3. Run int_defects_read_J_at_inc.
%
%% Check optional input argument
%Define tolerance on the contour integral path-dependence check.
if nargin == 2
    tol = 0.01; %Default of 1% difference between last and penultimate values
elseif nargin == 3
    tol = varargin{1};
else
    error('Unexpected number of input arguments.');
end

%% Read contour integral results from the .dat file
datScanInputStruct.filename = 'IntCrackJob1.dat';
[ outputArray, stepIncTimes ] = read_dat_contour_integral( datScanInputStruct, 'j' );

%% J-integral results
for k1 = 1:length(nodeNos)
    if isnan(nodeNos(k1))
        %If nodeNo is given as NaN, determine the greatest J-integral on the
        %crack tip line and return that.
        JLimsAll(1,:) = outputArray{2,maxInc-1}(end,:);
        JLimsAll(2,:) = outputArray{2,maxInc}(end,:);
        
        JAll = mean(JLimsAll,1);   %Determine node which has the highest J.
        [~,nodeNos(k1)] = max(JAll);    %Now assign nodeNos(k1) this value.
        
        JLims(k1,1) = JLimsAll(1,nodeNos(k1));     %J at minInc
        JLims(k1,2) = JLimsAll(2,nodeNos(k1));     %J at maxInc
        J(k1,1) = mean(JLims(k1,:));    %Best-guess J is simply the mean of the two.
    else
        %If nodeNo is given as a number, extract the J-integral for that node.
        JLims(k1,1) = outputArray{2,maxInc-1}(end,nodeNos(k1));     %J at minInc
        JLims(k1,2) = outputArray{2,maxInc}(end,nodeNos(k1));       %J at maxInc
        J(k1,1) = mean(JLims(k1,:));    %Best-guess J is simply the mean of the two.
    end
    
    %Model time at the moment of the best-guess J-integral
    modelTime = mean([stepIncTimes(maxInc-1,6),stepIncTimes(maxInc,6)]);
end

%% Contour-dependence check
for k1 = 1:length(nodeNos)
    pathDepWarnFlag(k1,1) = false;
    
    JMinVect = outputArray{2,maxInc-1}(:,nodeNos(k1));     %J at minInc
    JMaxVect = outputArray{2,maxInc}(:,nodeNos(k1));       %J at maxInc
    
    JMinRelDiff = abs((JMinVect(end)-JMinVect(end-1))./JMinVect(end));  %Relative difference between last and penultimate values
    JMaxRelDiff = abs((JMaxVect(end)-JMaxVect(end-1))./JMaxVect(end));
    
    if JMinRelDiff > tol || JMaxRelDiff > tol
        pathDepWarnFlag(k1,1) = true;
    end
end

end