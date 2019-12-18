function [] = int_defects_write_logfile_1( paramLogStruct )
% int_defects_write_logfile_1.m
% Harry Coules 2016
%
%DESCRIPTION
%This function writes a logfile containing information on the generation on
%Abaqus .inp files during studies of interacting defects. It is called at
%the end of int_defects_write_inp_parametric.m, and writes the logfile into
%the working folder.
%
%INPUT ARGUMENTS
% paramLogStruct - Structure containing a log of the parameters used for
%   the models.
%
%OUTPUT ARGUMENTS
%   *none*
%
%% Determine the numbers of failed .inp files and the reasons for them
%Scan paramLogStruct
for k1 = 1:size(paramLogStruct,1)
    inpStatus(k1) = paramLogStruct(k1).inpStatus;
    notValid(k1) = any(strfind(paramLogStruct(k1).inpWriteExitMessage,'modelParamStruct indicates this set of parameters not valid.'));
    shellSweepFailed(k1) = any(strfind(paramLogStruct(k1).inpWriteExitMessage,'Shell sweep feature failed'));
    featureCreationFailed(k1) = any(strfind(paramLogStruct(k1).inpWriteExitMessage,'Feature creation failed.'));
    notMeshed(k1) = any(strfind(paramLogStruct(k1).inpWriteExitMessage,'The following part instances are not meshed'));
    contourIntegralObjectProblem(k1) = any(strfind(paramLogStruct(k1).inpWriteExitMessage,'Warning: The following contour integral objects contains'));
    unknownReason(k1) = any(strfind(paramLogStruct(k1).inpWriteExitMessage,'The input file was not generated.'))|any(strfind(paramLogStruct(k1).inpWriteExitMessage,'#* Error:'));
end
noModels = length(inpStatus);
noInpStatus = sum(inpStatus);
noFailed = noModels-noInpStatus;

noNotValid = sum(notValid);
noShellSweepFailed = sum(shellSweepFailed);
noFeatureCreationFailed = sum(featureCreationFailed);
noNotMeshed = sum(notMeshed);
noContourIntegralObjectProblem = sum(contourIntegralObjectProblem);
noUnknownReason = sum(unknownReason);

%% Write the logfile
fid = fopen('log1.txt','at');
fprintf(fid, ' \n');
fprintf(fid, 'log1.txt\n');
fprintf(fid, [datestr(now),'\n']);
fprintf(fid, '***\n');
fprintf(fid, ['- inpStatus: ',num2str(noInpStatus),' .inp files written for ',num2str(noModels),' models.\n']);
fprintf(fid, '- Of the failures:\n');
fprintf(fid, [num2str(noNotValid),'/',num2str(noFailed),' - Model parameters not valid.\n']);
fprintf(fid, [num2str(noShellSweepFailed),'/',num2str(noFailed),' - Shell sweep feature creation failed.\n']);
fprintf(fid, [num2str(noFeatureCreationFailed),'/',num2str(noFailed),' - General feature creation failed.\n']);
fprintf(fid, [num2str(noNotMeshed),'/',num2str(noFailed),' - Meshing failed.\n']);
fprintf(fid, [num2str(noContourIntegralObjectProblem),'/',num2str(noFailed),' - Problem with creation of contour integral objects.\n']);
fprintf(fid, [num2str(noUnknownReason),'/',num2str(noFailed),' - Unknown cause of failure.\n']);
fprintf(fid, '***\n');
fclose all;

end

