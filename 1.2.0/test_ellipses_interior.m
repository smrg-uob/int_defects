function [ interiorFlag ] = test_ellipses_interior( xEC, yEC, r1, r2, xTP, yTP )
%test_ellipses_interior.m
%Harry Coules 2017
%
%DESCRIPTION
%Simple test for whether a point lies within an ellipse.
%
%INPUT ARGUMENTS
%   xEC = x-coord. of the centre of the ellipse.
%   yEC = y-coord. of the centre of the ellipse.
%   r1 = radius of the ellipse in x.
%   r2 = radius of the ellipse in y.
%   xTP = x-coord. of test point.
%   yTP = y-coord. of test point.
%
%OUTPUT ARGUMENTS
%   interiorFlag = Logical true if test point is inside the ellipse, false
%       otherwise.

%%
a = (((xTP-xEC).^2)./(r1.^2))+(((yTP-yEC).^2)./(r2.^2));
if a <= 1
    interiorFlag = true;
else
    interiorFlag = false;
end

end

