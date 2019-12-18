function [varargout] = plot_KI_over_crack( paramLogStruct, outputArray, modelNo, crackNo, contourNo, varargin )
%plot_KI_over_crack.m
%Harry Coules 2015
%
%DESCRIPTION
%This function generates a plot of KI vs ellipse parametric angle for an
%elliptical or semi-elliptical crack-like defect in a structure. It uses
%data from the input data structure modelParamStruct, which contains
%results from a set of finite elements created and run using the
%int_defects package.
%
%INPUT ARGUMENTS
%   modelParamStruct - Structure containing model parameters
%   outputArray - Cell array containing results of the analysis
%   modelNo - Number of the model for which results are to be plotted
%   crackNo - Number of the crack for which results are to be plotted (1 for
%       Crack A, 2 for Crack B).
%   contourNo - Contour number to use
%
%OPTIONAL INPUT ARGUMENTS
%   lineColorStr - A string to define the colour of the plotted line.
%
%OPTIONAL OUTPUT ARGUMENT
%   lineHandle - The handle of the plotted line.
%
%% Line colour, if present
if isempty(varargin)
    lineColorStr = 'k';     %Default colour is black
elseif length(varargin)==1
    if isstr(varargin{1})
        lineColorStr = varargin{1};
    else
        error('If line colour is specified as an input argument, it must be a string.');
    end
else
    error('Unexpected number of input arguments');
end

%% Get data to plot
%Extract SIF and arc length from array
k1 = outputArray{modelNo}.k{2,crackNo}(contourNo,:);    %Mode I SIF
lOverL = (0:(length(k1)-1))/(length(k1)-1);                         %Normalised arc length

%Crack sizes
if crackNo==1
    a1 = paramLogStruct(modelNo).modelParams.crackSizes.aA1;
    a2 = paramLogStruct(modelNo).modelParams.crackSizes.aA2;
    
    %For Crack A, check for flag indicating that it's a subsurface crack.
    if isfield(paramLogStruct(modelNo).modelParams.crackPositions,'subsurfaceAFlag')
        semiEllipseFlag = ~paramLogStruct(modelNo).modelParams.crackPositions.subsurfaceAFlag;
    else
        semiEllipseFlag = true; %Assume crack is on the surface unless otherwise specified.
    end
elseif crackNo ==2
    a1 = paramLogStruct(modelNo).modelParams.crackSizes.aB1;
    a2 = paramLogStruct(modelNo).modelParams.crackSizes.aB2;
    
    warning('Normalised arc length is flipped for Crack B.')
    lOverL = fliplr(lOverL);    %lOverL needs to be flipped due to different ordering of crack tip nodes for Crack B
    
    %For Crack B, check for flag indicating that it's a subsurface crack.
    if isfield(paramLogStruct(modelNo).modelParams.crackPositions,'subsurfaceBFlag')
        semiEllipseFlag = ~paramLogStruct(modelNo).modelParams.crackPositions.subsurfaceBFlag;
    else
        semiEllipseFlag = true; %Assume crack is on the surface unless otherwise specified.
    end
else
    error('Invalid input for crackNo. Should be 1 (for Crack A) or 2 (for Crack B).')
end

%Calculate ellipse parametric angles
[ t ] = calc_ellipse_angle( a2, a1, lOverL, true, semiEllipseFlag );
tDeg = 180*t/pi;

%% Plotting
lineHandle = plot(tDeg, k1,'color',lineColorStr,'linewidth',1);
hold on;

%% Optional output
if nargout == 1
    varargout{1} = lineHandle;
elseif nargout ~= 0
    error('Unexpected number of output arguments.');
end

end