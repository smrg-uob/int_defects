function [] = int_defects_plot_J_vs_stress( filename, nodeNo )
%int_defects_plot_J_vs_stress.m
%Harry Coules 2018
%
%DESCRIPTION
%This function determines the J-integral and stress for a given position
%on the crack tip line for a cracked-body FE model. It then plots J vs
%stress. NOTE: This function assumes that the stress ramps linearly with
%the model step time.
%
%INPUT ARGUMENTS
%   filename - Name of the .dat file from which results should be taken.
%   nodeNo - Node number specifying the location on the crack tip line
%       which is being interrogated.
%
%OUTPUT ARGUMENTS
%   J - J-integral at the requested node no.
%   s - Corresponding remote stress.
%
%% Read data from the .dat file
[ outputArray, stepIncTimes ] = read_dat_contour_integral( filename, 'j' );

%% Determine J-integral at the desired point, and the stress
s = stepIncTimes(:,6);  %Note that this assumes that the stress is ramped linearly with step time.

J = [];
for k1 = 1:size(stepIncTimes,1)
    J(k1) = outputArray{2,k1}(end,nodeNo);
end

%% Plot
hold on
plot(s,J,'.-');

end