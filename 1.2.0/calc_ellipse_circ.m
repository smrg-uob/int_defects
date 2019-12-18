function [ circ ] = calc_ellipse_circ(a, b, varargin)
%calc_ellipse_circ.m
%Harry Coules 2019
%
%DESCRIPTION
%This function approximates the circumference of an ellipse from the
%lengths of its semi-major axes. It uses the (exact) infinite sum formula
%determined by James Ivory and Friedrich Bessel.
%
%INPUT ARGUMENTS
%   a - Length of first semi-major axis.
%   b - Length of second semi-major axis.
% *OPTIONAL*
%   n - Number of terms of the infinite sum to use (default = 100 if
%       input argument is omitted).
%
%OUTPUT ARGUMENTS
%   circ - Circumference of the ellipse.
%
%% Optional input argument
if nargin==2
    n = 100;    %Silently default to n=100
elseif nargin==3
    n = varargin{1};
elseif nargin>3
    error('Too many input arguments.');
end

%% Calculate circumference
h = ((a-b).^2.)/((a+b).^2);

circ = pi*(a+b);    %0th order coefficient
for k1 = 1:n
    binCoeff = gamma(1.5)./(gamma(1+k1).*gamma(1.5-k1));  %Gamma function formula for binomial coefficient
    circ = circ + pi*(a+b).*(h.^k1).*(binCoeff.^2);       %Use binoomial coeffient in Ivory/Bessel formula
end

end
