function [ intersectionFlag, phi, sA, sB, sMin, indexA, indexB ] = calc_ellipses_approach( rA1, rA2, rB1, rB2, d, b )
%calc_ellipses_approach.m
%Harry Coules 2015
%
%DESCRIPTION
%This function approximates the distance between two ellipses. It is
%intended to be used for studying the interaction between a pair of
%semi-elliptical defects (or an semi-elliptical surface defect and a buried)
%elliptical defect) in a plate.
% - 12/12/17: Modified to use complete ellipses instead of semi-ellipses,
%   since this was causing problems for buried elliptical cracks. Also
%   added function to check whether the centre of each ellipse lies within
%   the other.
%
%INPUT ARGUMENTS
%   rA1 -  Radius of Ellipse A in the direction of semi-axis 1 (i.e. the
%       through-thickness or y direction).
%   rA2 - Radius of Ellipse A in the direction of semi-axis 2 (i.e. x,
%       along the plate's surface).
%   rB1 - Radius 1 of Ellipse B
%   rB2 - Radius 2 of Ellipse B
%   d - Distance in x between the ends of each ellipse. If this is negative
%       then the end of each ellipse overhangs the other one.
%   b - Distance in y of the centres of each semi-ellipse from each other.
%       For semi-elliptical cracks emanating from opposite sides of a
%       plate, this is equal to the thickness of the plate.
%
%OUTPUT ARGUMENTS
%   intersectionFlag - Logical set to true if the ellipses intersect (or
%       if one fully encloses the another).
%   phi - Set of ellipse parametric angles used for numerical approximation
%       of approach distances.
%   sA - Vector of distances to Ellipse B for each point in phi around
%       Ellipse A.
%   sB - Vector of distances to Ellipse A for each point in phi around
%       Ellipse B.
%   sMin - Closest approach between the ellipses.
%   indexA - Index for the phi angle on Ellipse A at which closest approach
%       occurs.
%   indexB - Index for the phi angle on Ellipse B at which closest approach
%       occurs.
%
%% Calculate coordinates of the lines of the two ellipses
%Ellipse parametric angle
phi = linspace(0,2*pi,2001)';

%Calculate coordinates for Ellipse A
xA = rA2*cos(phi);
yA = rA1*sin(phi);

%Calculate coordinates for Ellipse B - note centre is at [cB, b]
cB = rA2 + rB2 + d; %Centre of Ellipse B
xB = rB2*cos(phi);
xB = cB - xB;
yB = b - rB1*sin(phi);

%% Check whether the centre of each ellipse is enclosed by the other ellipse
[ interiorFlag1 ] = test_ellipses_interior( 0, 0, rA2, rA1, cB, b );  %Check if ellipse B centre is inside ellipse A.
[ interiorFlag2 ] = test_ellipses_interior( cB, b, rB2, rB1, 0, 0 );  %Check if ellipse A centre is inside ellipse B.
interiorFlag = interiorFlag1|interiorFlag2;

%% Check for intersection of the ellipse lines
%Note that intersectionFlag is set as true when the centre of one ellipse
%is inside the other, regardless of whether the bounding lines ever cross.
%This could occur if one ellipse is completely enclosed by the other.
[xInt,~] = intersections(xA,yA,xB,yB);
if isempty(xInt) && ~interiorFlag
    intersectionFlag = false;
else
    intersectionFlag = true;
end

%% If intersection does not occur, calculate closest approaches
s = [];
sA = [];
sB = [];
sMin = NaN;
indexA = NaN;
indexB = NaN;
if ~intersectionFlag
    for k1 = 1:length(xA);
        s(k1,:) = ((((xB-xA(k1)).^2)+((yB-yA(k1)).^2)).^0.5)';
    end
    
    %Calculate approach distances
    sA = min(s,[],2);  %The distance to Ellipse B for each point on Ellipse A
    sB = min(s,[],1);  %The distance to Ellipse A for each point on Ellipse B
    sB = sB';
    
    %Calculate closest approach and the index at which it occurs
    [sAmin, indexA] = min(sA);
    [sBmin, indexB] = min(sB);
    
    %Check that the minimum distance between ellipses is the same when
    %calculated starting from each ellipse.
    if sAmin ~= sBmin
        warning('Closest approach distances calculated from each of the semi-ellipses are not the same.')
    end
    sMin = sAmin;
end

end