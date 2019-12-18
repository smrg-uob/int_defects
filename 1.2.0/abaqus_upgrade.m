function [ ] = abaqus_upgrade( jobName )
%abaqus_upgrade.m
%Harry Coules 2018
%
%DESCRIPTION
%This function upgrades an Abaqus .odb file to be compatible with the
%latest version of Abaqus/CAE, using Abaqus' input file and output database
%upgrade utility.
%
%INPUT ARGUMENTS
%   jobName - Name of the Abaqus job for which there is an odb file to be
%       upgraded.
%
%OUTPUT ARGUMENTS
%   *None*
%   (.odb file is modified)
%
%%
%Check if we're running on a PC or a unix system
if ispc
    PCflag = true;
elseif isunix
    PCflag = false;
else
    error('Unrecognised system type - abaqus_run only works on PCs and unix systems.')
end

%Run upgrade
commandStr = ['abaqus upgrade job=',jobName,'_UPGR odb=',jobName,' interactive'];
if PCflag
    [~,~] = dos(commandStr,'-echo');
elseif ~PCflag
    [~,~] = unix(commandStr,'-echo');
else
    error('Unrecognised system type - abaqus_run only works on PCs and unix systems.')
end

%Delete old .odb file and rename the new one
delete([jobName,'.odb']);
pause(5);
movefile([jobName,'_UPGR.odb'],[jobName,'.odb']);

end
