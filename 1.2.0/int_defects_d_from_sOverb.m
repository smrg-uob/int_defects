function [ d, varargout] = int_defects_d_from_sOverb( aA1, aA2, aB1, aB2, b, sOverb, varargin )
%int_defects_d_from_sOverb.m
%Harry Coules 2016
%
%DESCRIPTION
%This function is used to approximate the distance d separating the
%edges of two through-thickness semi-elliptical cracks emanating from
%opposite sides of a finite-thickness plate, from the (normalised) closest
%approach of the crack fronts sOverb.
%
%INPUT ARGUMENTS
%   aA1 - Depth of crack A.
%   aA2 - Half-width of crack A.
%   aB1 - Depth of crack B.
%   aB2 - Half-width of crack B.
%   b - Thickness of the plate.
%   sOverb - Closest approach of crack fronts.
%OPTIONAL INPUT ARGUMENT
%   tol - Relative tolerance on the stopping criterion for the iterative
%         method. It stops when it is within tol*sOverb of the required value
%         of sOverb. If omitted, the default is 1e-3.
%
%OUTPUT ARGUMENTS
%   d - The distance between crack edges in the along-plate direction. If
%       the required sOverb cannot be achieved, d is given as NaN and a
%       warning is shown.
%   *OPTIONAL*
%   noIters - The number of iterations performed to determine d.
%
%% Check optional input argument
if isempty(varargin)
    tol = 1e-3;
elseif length(varargin) == 1
    tol = varargin{1};
else
    error('Too many input arguments.')
end

%% Iterative method to find d
%Determine the value of d required to achieve the desired sOverb
[ ~, ~, ~, ~, sMin, ~, ~ ] = calc_ellipses_approach( aA1, aA2, aB1, aB2, -aA2-aB2, b );  %First run calc_ellipses_approach for cracks directly opposite each other
if sMin > sOverb*b
    d = NaN;
    warning('Required separation is too small. It cannot be achieved even at zero offset for this crack pair. Returning d = NaN.')
else
    %Now iterate to find the correct d
    d = sOverb*b;  %First separation distance to try
    stepSize = sOverb*b; %Initial step size
    separationDone = false;
    noIters = 0;
    while ~separationDone
        noIters = noIters+1;
        [ intersectionFlag, ~, ~, ~, sMin, ~, ~ ] = calc_ellipses_approach( aA1, aA2, aB1, aB2, d, b );
        if intersectionFlag || sMin < sOverb*b || d < -aA2-aB2 %If we've gone too far...
            d = dOld;   %Go back to the last d - note this condition should never occur on the 1st iteration when there is no dOld.
            stepSize = 0.5*stepSize; %Halve the step size
        elseif (1-tol)*sOverb < abs(sMin/b) && (1+tol)*sOverb > abs(sMin/b) %Stopping criterion: we must be within tol*sOverb of correct value of sOverb
            separationDone = true;                                          %Note that we check for the stopping criterion AFTER checking for overshoot etc.
        else
            dOld = d;
            d = d - stepSize; %Otherwise modify d
        end
    end
end

%% Optional output argument
if nargout > 2
    error('Too many output arguments specified.')
elseif nargout == 2
    varargout{1} = noIters; %Number of iterations performed.
end

end