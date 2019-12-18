function [varargout] = plot_IF_over_crack( paramLogStruct, outputArray, modelNo, crackNo, contourNo, semiEllipseFlag, varargin )
%plot_IF_over_crack.m
%Harry Coules 2015
%
%This function generates a plot of interaction factor vs ellipse parametric
%angle for a single model of a pair of interacting elliptical or semi-
%elliptical crack-like defects in structure. It uses data from the data
%array outputArray, which contains results from a set of finite
%elements created and run using the int_defects package. outputArray is
%created by the function int_defects_calc_interaction_factor.m.
%
%INPUT ARGUMENTS
%   modelParamStruct - Structure containing model parameters
%   outputArray - Cell array containing results of the analysis
%   modelNo - Number of the model for which results are to be plotted
%   crackNo - Number of the crack for which results are to be plotted (1 for
%       Crack A, 2 for Crack B).
%   contourNo - Contour number to use.
%   semiEllipseFlag - Logical indicating whether the crack is a
%       semi-ellipse or a full ellipse (i.e. a subsurface defect).
%
%OPTIONAL INPUT ARGUMENTS
%   lineColorStr - A string to define the colour of the plotted line.
%   flipPhiFlag - A logical indicating whether phi chould be flipped for
%   plotting.
%
%OPTIONAL OUTPUT ARGUMENT
%   lineHandle - The handle of the plotted line.
%
%% Line colour, if present
if isempty(varargin)
    lineColorStr = 'k';     %Default colour is black
    flipPhiFlag = false;    %By default do not flip phi
elseif length(varargin)==1
    if isstr(varargin{1})
        lineColorStr = varargin{1};
        flipPhiFlag = false;
    elseif islogical(varargin{1})
        flipPhiFlag = varargin{1};
        lineColorStr = 'k';
    else
        error('If line colour is specified as an input argument, it must be a string.');
    end
elseif length(varargin)==2
    if isstr(varargin{1}) && islogical(varargin{2})
        lineColorStr = varargin{1};
        flipPhiFlag = varargin{2};
    else
        error('If both optional input arguments are specified, the first must be a string and the second must be a logical.')
    end
else
    error('Unexpected number of input arguments');
end

%% Get data to plot
%Extract SIF and arc length from array
intFactor = outputArray{modelNo}.interactionFactor{1,crackNo}(contourNo,:);    %Interaction factor
lOverL = (0:(length(intFactor)-1))/(length(intFactor)-1);                      %Normalised arc length

%Crack sizes
if crackNo==1
    a1 = paramLogStruct(modelNo).modelParams.crackSizes.aA1;
    a2 = paramLogStruct(modelNo).modelParams.crackSizes.aA2;
elseif crackNo ==2
    a1 = paramLogStruct(modelNo).modelParams.crackSizes.aB1;
    a2 = paramLogStruct(modelNo).modelParams.crackSizes.aB2;
    
    warning('Normalised arc length is flipped for Crack B.')
    lOverL = fliplr(lOverL);    %lOverL needs to be flipped due to different ordering of crack tip nodes for Crack B
else
    error('Invalid input for crackNo. Should be 1 (for Crack A) or 2 (for Crack B).')
end

%Calculate ellipse parametric angles
[ phi ] = calc_ellipse_angle( a2, a1, lOverL, true, semiEllipseFlag, 1e6);
if flipPhiFlag
    phi = flipud(phi);
end
phiDeg = 180*phi/pi;
phiNorm = phi/pi;

%% Plotting
lineHandle = plot(phiDeg, intFactor,'color',lineColorStr,'linewidth',1);
hold on;

%% Optional output
if nargout == 1
    varargout{1} = lineHandle;
elseif nargout ~= 0
    error('Unexpected number of output arguments.');
end

end