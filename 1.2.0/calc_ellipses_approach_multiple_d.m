function [ intersectionFlag, sMin ] = calc_ellipses_approach_multiple_d( rA1, rA2, rB1, rB2, d, b )
%calc_ellipses_approach_multiple_d.m
%Harry Coules 2015
%
%This function is used to approximate the distance between two semi-
%ellipses spaced at a range of distances. It runs calc_ellipses_approach
%multiple times, once for each separation distance, and reports only
%whether the ellipses intersect and if not, the minimum distance between
%them.
%
%INPUT ARGUIMENTS
%   rA1 -  Radius of Ellipse A in the direction of semi-axis 1 (i.e. the
%       through-thickness or y direction).
%   rA2 - Radius of Ellipse A in the direction of semi-axis 2 (i.e. x,
%       along the plate's surface).
%   rB1 - Radius 1 of Ellipse B
%   rB2 - Radius 2 of Ellipse B
%   d - Distance(s) in x between the ends of each ellipse. If this is negative
%       then the end of each ellipse overhangs the other one. NOTE: d can
%       be a vector.
%   b - Distance in y of the centres of each semi-ellipse from each other.
%       For semi-elliptical cracks emanating from opposite sides of a
%       plate, this is equal to the thickness of the plate.
%
%OUTPUT ARGUMENTS
%   intersectionFlag - Logical set to true if the ellipses intersect, or if
%       one fully encloses the another. This is a vector if d is a vector.
%   sMin - Closest approach between the ellipses (as a vector if d is a
%       vector).
%
%%
for k1 = 1:length(d)
    [ intersectionFlag(k1),~,~,~,sMin(k1),~,~] = calc_ellipses_approach( rA1, rA2, rB1, rB2, d(k1), b );
end

end