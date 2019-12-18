function [] = int_defects_write_logfile_2( mainStruct, singleStruct, varargin)
%int_defects_write_logfile_2.m
%Harry Coules 2016
%
%DESCRIPTION
%This function writes a logfile containing information about the status of
%the models in a parametric series, and whether the results from these
%models were successfully used to calculate the interaction factors for
%pairs of cracks.
%
%INPUT ARGUMENTS
% mainStruct - Structure containing information about the main set of
%   models, i.e. the crack-pair models.
% singleStruct - Structure containing information about the (first) set of
%   single crack models.
% *OPTIONAL INPUT ARGUMENT*
% singleStruct2 - Structure containing information about the second set of
%   single crack models, if a second set was used.
%
%OUTPUT ARGUMENTS
% *none*

%% Optional input argument
if length(varargin) == 1
    singleStruct2 = varargin{1};
elseif length(varargin) > 1
    error('No more than three input arguments should be provided to this function.');
end

%% Determine the numbers of failed models the reasons for them
%Determine the number of parameter sets where interaction factors couldn't
%be calculated.
for k1 = 1:size(mainStruct.outputArray,2)
    intFactorsCalced(k1) = mainStruct.outputArray{1,k1}.intFactorsCalced;
    singleIntFactorsCalced(k1,:) = mainStruct.outputArray{1,k1}.singleIntFactorsCalced;
end

totalModels = size(mainStruct.outputArray,2);
noIntFactorsCalced = sum(intFactorsCalced);
noIntFactorsFailed = totalModels - noIntFactorsCalced;

%Determine the cause of failed interaction factor calcs
    mainErrors = 0;
    AErrors = 0;
    BErrors = 0;
    multipleErrors = 0;
    
    mainNoInp = 0;
    mainJobFailed = 0;
    mainUnknownFail = 0;
    
    ANoInp = 0;
    AJobFailed = 0;
    AUnknownFail = 0;
    
    BNoInp = 0;
    BJobFailed = 0;
    BUnknownFail = 0;
    
    multipleNoInp = 0;
    multipleJobFailed = 0;
    multipleUnknownFail = 0;
    
    IFCalcErrors = 0;
    
    statusMain = [];
    statusA = [];
    statusB = [];
for k1 = 1:size(mainStruct.paramLogStruct,1)
    %Single crack model numbers
    singleCrackNoA(k1) = mainStruct.paramLogStruct(k1).singleCrackNo(1);
    singleCrackNoB(k1) = mainStruct.paramLogStruct(k1).singleCrackNo(2);
    
    %Status of main model
    statusMain(k1) = mainStruct.outputArray{k1}.datReadSuccessFlag;
    exitStatusMain(k1) = mainStruct.paramLogStruct(k1).jobExitStatus;
    
    %Status of Model A
    statusA(k1) = singleStruct.outputArray{singleCrackNoA(k1)}.datReadSuccessFlag;
    exitStatusA(k1) = singleStruct.paramLogStruct(singleCrackNoA(k1)).jobExitStatus;
    
    %Status of Model B
    if exist('singleStruct2','var')
        statusB(k1) = singleStruct2.outputArray{singleCrackNoB(k1)}.datReadSuccessFlag;
        exitStatusB(k1) = singleStruct2.paramLogStruct(singleCrackNoB(k1)).jobExitStatus;
    else
        statusB(k1) = singleStruct.outputArray{singleCrackNoB(k1)}.datReadSuccessFlag;
        exitStatusB(k1) = singleStruct.paramLogStruct(singleCrackNoB(k1)).jobExitStatus;
    end
    
    %Check reason
    if statusMain(k1)==1 && statusA(k1)==1 && statusB(k1)==1
        if ~intFactorsCalced(k1)
            IFCalcErrors = IFCalcErrors+1;
        end
    elseif ~statusMain(k1) && statusA(k1) && statusB(k1)
        if intFactorsCalced(k1)
            error(['Interaction factors appear to have been calculated for a case for which not all model results are available. k1 = ',num2str(k1)]);
        end
        
        mainErrors = mainErrors+1;
        
        if exitStatusMain(k1) == 0
            mainJobFailed = mainJobFailed+1;
        elseif exitStatusMain(k1) == 2
            mainNoInp = mainNoInp+1;
        else
            mainUnknownFail = mainUnknownFail+1;
        end
    elseif statusMain(k1) && ~statusA(k1) && statusB(k1)
        if intFactorsCalced(k1)
            error(['Interaction factors appear to have been calculated for a case for which not all model results are available. k1 = ',num2str(k1)]);
        end       
        
        AErrors = AErrors+1;
        
        if exitStatusA(k1) == 0
            AJobFailed = AJobFailed+1;
        elseif exitStatusA(k1) == 2
            ANoInp = ANoInp+1;
        else
            AUnknownFail = AUnknownFail+1;
        end
    elseif statusMain(k1) && statusA(k1) && ~statusB(k1)
        if intFactorsCalced(k1)
            error(['Interaction factors appear to have been calculated for a case for which not all model results are available. k1 = ',num2str(k1)]);
        end
        
        BErrors = BErrors+1;
        
        if exitStatusB(k1) == 0
            BJobFailed = BJobFailed+1;
        elseif exitStatusB(k1) == 2
            BNoInp = BNoInp+1;
        else
            BUnknownFail = BUnknownFail+1;
        end
    else
        if intFactorsCalced(k1)
            error(['Interaction factors appear to have been calculated for a case for which not all model results are available. k1 = ',num2str(k1)]);
        end
        
        multipleErrors = multipleErrors+1;
        
        if (exitStatusMain(k1)==0)+(exitStatusA(k1)==0)+(exitStatusB(k1)==0)>1
            multipleJobFailed = multipleJobFailed+1;
        elseif (exitStatusMain(k1)==2)+(exitStatusA(k1)==2)+(exitStatusB(k1)==2)>1
            multipleNoInp = multipleNoInp+1;
        else
            multipleUnknownFail = multipleUnknownFail+1;
        end
    end 
end

%% Write the logfile
fid = fopen('log2.txt','at');
fprintf(fid, ' \n');
fprintf(fid, 'log2.txt\n');
fprintf(fid, [datestr(now),'\n']);
fprintf(fid, '***\n');

%No. of failures
fprintf(fid, ['- Full IF data calculated for ',num2str(noIntFactorsCalced),' out of ',num2str(totalModels),' parameter sets.\n']);
fprintf(fid, '- Of the failures:\n');

%IF calculation error
fprintf(fid, [num2str(IFCalcErrors),'/',num2str(noIntFactorsFailed),' - All relevant models ran, but unable to calculate IFs.\n']);

%Multiple errors
fprintf(fid, [num2str(multipleErrors),'/',num2str(noIntFactorsFailed),' - Multiple model results not available, of which:\n']);
fprintf(fid, ['   ',num2str(multipleNoInp),'/',num2str(multipleErrors),' - Multiple models (A/B/interacting) missing valid .inp files (model parameters probably invalid).\n']);
fprintf(fid, ['   ',num2str(multipleJobFailed),'/',num2str(multipleErrors),' - Multiple models (A/B/interacting) attempted but failed.\n']);
fprintf(fid, ['   ',num2str(multipleUnknownFail),'/',num2str(multipleErrors),' - Differing and/or unknown causes of failure.\n']);

%Main model
fprintf(fid, [num2str(mainErrors),'/',num2str(noIntFactorsFailed),' - Interacting crack model result not available, of which:\n']);
fprintf(fid, ['   ',num2str(mainNoInp),'/',num2str(mainErrors),' - No valid .inp file available for interacting crack model.\n']);
fprintf(fid, ['   ',num2str(mainJobFailed),'/',num2str(mainErrors),' - Interacting crack model attempted but failed.\n']);
fprintf(fid, ['   ',num2str(mainUnknownFail),'/',num2str(mainErrors),' - Unknown cause of failure.\n']);

%Crack A model
fprintf(fid, [num2str(AErrors),'/',num2str(noIntFactorsFailed),' - Crack A model result not available, of which:\n']);
fprintf(fid, ['   ',num2str(ANoInp),'/',num2str(AErrors),' - No valid .inp file available for Crack A model.\n']);
fprintf(fid, ['   ',num2str(AJobFailed),'/',num2str(AErrors),' - Crack A model attempted but failed.\n']);
fprintf(fid, ['   ',num2str(AUnknownFail),'/',num2str(AErrors),' - Unknown cause of failure.\n']);

%Crack B model
fprintf(fid, [num2str(BErrors),'/',num2str(noIntFactorsFailed),' - Crack B model result not available.\n']);
fprintf(fid, ['   ',num2str(BNoInp),'/',num2str(BErrors),' - No valid .inp file available for Crack B model.\n']);
fprintf(fid, ['   ',num2str(BJobFailed),'/',num2str(BErrors),' - Crack B model attempted but failed.\n']);
fprintf(fid, ['   ',num2str(BUnknownFail),'/',num2str(BErrors),' - Unknown cause of failure.\n']);

fprintf(fid, '***\n');
fclose all;

end

