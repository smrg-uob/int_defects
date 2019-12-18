function [ ] = int_defects_upgrade_odbs(jobName)
%int_defects_upgrade_odbs.m
%Harry Coules 2018
%
%DESCRIPTION
%This function loops through a set of numerically-named directories
%which each contain the results from a single Abaqus model, and upgrades
%any .odb files therein to be compatible with the latest installed version
%of Abaqus/CAE.
%
%INPUT ARGUMENTS
%   jobName - Name of the Abaqus jobs to be upgraded. Must be the same for
%       all directories.
%
%OUTPUT ARGUMENTS
% *none*
%
%USAGE
% 1. Navigate to the directory containing the sequentially-numbered
%   subdirectories which contain .odb files.
% 2. Call this function.
%%
k1 = 1;
while logical(exist(num2str(k1,'%06i'),'dir'))  %Note that use of this this while loop assumes that the directories are sequentially named with no ommissions.
    disp(['Upgrading .odb file from directory ',num2str(k1,'%06i')]);
    tic;
    cd(num2str(k1,'%06i'));
    abaqus_upgrade(jobName);   %Run abaqus_cleanup targetting .odb files only
    cd ..
    disp(['Done. Time taken for this .odb file: ',num2str(toc),'seconds.']);
    k1 = k1+1;
end

end
