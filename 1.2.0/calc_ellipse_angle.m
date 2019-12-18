function [ t ] = calc_ellipse_angle( a1, a2, lOverL, varargin )
%calc_ellipse_angle.m
%Harry Coules 2015
%
%This function is used to approximate the ellipse parametric angle
%corresponding to a particular length of arc around the ellipse's
%circumference, measured from the end of one of its axes.
%
%INPUT ARGUMENTS
%   a1 - Ellipse axis 1 length (axis from which t is measured)
%   a2 - Ellipse axis 2 length
%   lOverL - Normalised path length around the perimeter of the ellipse
%   (normalised to the path length around one half of the ellipse). If the
%   optional input argument lengthIsNormalised is 'false' then this can be
%   provided in un-normalised form. This is measured from Ellipse axis 1.
% OPTIONAL INPUT ARGUMENTS
%   lengthIsNormalised - Logical to indicate whether the third input
%       argument is normalised. Default is true.
%   semiEllipseFlag - Logical to indicate whether the geometry is a
%       semi-ellipse rather than a complete ellipse. Default is true.
%   noChords - Number of chords to use in in the numerical approximation of
%       the angle. Default is 1e6.
%
%OUTPUT ARGUMENTS
%   t - Ellipse parametric angle
%
%NOTES
% - This function supersedes approx_ellipse_angle.m.
% - When optional input arguments are given, they must be in the order
%   listed above. When semiEllipseFlag is given it must be preceeded by
%   lengthIsNormalised. When noChords is given it must be preceeded by
%   lengthIsNormalised and semiEllipseFlag.
%
%% Optional input arguments
if isempty(varargin)
    %Default to assuming that the length is normalised to the total length 
    %of a semi-ellipse. Note that there is no warning for this.
    lengthIsNormalised = true;
    semiEllipseFlag = true;
    noChords = 1e6;
elseif length(varargin) == 1
    lengthIsNormalised = varargin{1};
    semiEllipseFlag = true;
    noChords = 1e6;
elseif length(varargin) == 2
    lengthIsNormalised = varargin{1};
    semiEllipseFlag = varargin{2};
    noChords = 1e6;
elseif length(varargin) == 2
    lengthIsNormalised = varargin{1};
    semiEllipseFlag = varargin{2};
    noChords = 1e6;
elseif length(varargin) == 3
    lengthIsNormalised = varargin{1};
    semiEllipseFlag = varargin{2};
    noChords = varargin{3};
else
    error('Too many input arguments.')
end

%% Numerical integration to get an array of arc lengths for different parametric angles
%Vector of ellipse parametric angles
if semiEllipseFlag
    ts = linspace(0,pi,noChords);
else
    ts = linspace(0,2*pi,noChords);
end

%Length of each chord
integrand = sqrt(((a1.^2).*(sin(ts).^2))+...
    ((a2.^2).*(cos(ts).^2)));

%Integrate
ls = cumtrapz(ts,integrand);

%% Use total length of the semi-ellipse to get absolute path length
L = ls(end);
if lengthIsNormalised
    l = lOverL.*L;  %If the user has provided the normalised path length as input
else
    l = lOverL; %If the user has provided the absolute path length as input
end

%% Find the parametric angle(s) with the arc length(s) closest to what we want
t = zeros(length(l),1);
for k1 = 1:length(l)
    [~,index] = min(abs(ls-l(k1)));
    t(k1) = ts(index);
end

end