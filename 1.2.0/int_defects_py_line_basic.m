function [ lineOut ] = int_defects_py_line_basic( lineIn, modelParamStruct )
%int_defects_py_line_basic.m
%Harry Coules 2015
%
%DESCRIPTION
%This function is used as part of the line-by-line writing of a python
%script, which is then used by the Abaqus/CAE scripting interface to
%generate a model. It is part of the int_defects software package for
%analysing structures containing interacting crack-like defects.
%
%This function scans the string in lineIn and modifies it based on
%its contents, adding data defined in modelParamStruct. It is called by
%int_defects_write.py.
%
%INPUT ARGUMENTS
%   lineIn - String containing the current line read from the master .py
%       file.
%   modelParamStruct - Structure containing model data, which will be used
%       to modify lineIn.
%
%OUTPUT ARGUMENTS
%   lineOut - String containing a modified version of lineIn, with
%       modelling parameters from modelParamStruct added, if necessary.
%
%NOTES
%   - The function int_defects_py_line_edges.m is also used to make edits
%       to the .py file, in a similar manner to this function.
%
%%
line = lineIn;
formatSpec = '%.7f';   %Specify the number format to be written.

%%
%Crack Profile sizes
if any(strfind(line,'#aA1#'))
    paramInd = strfind(line,'#aA1#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.crackSizes.aA1,formatSpec);
    linePart3 = line((paramInd+5):end);
    line = [linePart1,linePart2,linePart3];
end

if any(strfind(line,'#aA2#'))
    paramInd = strfind(line,'#aA2#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.crackSizes.aA2,formatSpec);
    linePart3 = line((paramInd+5):end);
    line = [linePart1,linePart2,linePart3];
end

if any(strfind(line,'#aB1#'))
    paramInd = strfind(line,'#aB1#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.crackSizes.aB1,formatSpec);
    linePart3 = line((paramInd+5):end);
    line = [linePart1,linePart2,linePart3];
end

if any(strfind(line,'#aB2#'))
    paramInd = strfind(line,'#aB2#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.crackSizes.aB2,formatSpec);
    linePart3 = line((paramInd+5):end);
    line = [linePart1,linePart2,linePart3];
end

%Crack positions
if any(strfind(line,'#xA#'))
    paramInd = strfind(line,'#xA#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.crackPositions.xA,formatSpec);
    linePart3 = line((paramInd+4):end);
    line = [linePart1,linePart2,linePart3];
end

if any(strfind(line,'#xB#'))
    paramInd = strfind(line,'#xB#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.crackPositions.xB,formatSpec);
    linePart3 = line((paramInd+4):end);
    line = [linePart1,linePart2,linePart3];
end

if any(strfind(line,'#yA#'))
    paramInd = strfind(line,'#yA#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.crackPositions.yA,formatSpec);
    linePart3 = line((paramInd+4):end);
    line = [linePart1,linePart2,linePart3];
end

if any(strfind(line,'#yB#'))
    paramInd = strfind(line,'#yB#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.crackPositions.yB,formatSpec);
    linePart3 = line((paramInd+4):end);
    line = [linePart1,linePart2,linePart3];
end

if any(strfind(line,'#yBminusb#'))
    paramInd = strfind(line,'#yBminusb#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.crackPositions.yB-modelParamStruct.plateSizes.b,formatSpec);
    linePart3 = line((paramInd+10):end);
    line = [linePart1,linePart2,linePart3];
end

%Crack tip zone sizes
if any(strfind(line,'#rpA1#'))
    paramInd = strfind(line,'#rpA1#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.crackTipZoneSizes.rpA1,formatSpec);
    linePart3 = line((paramInd+6):end);
    line = [linePart1,linePart2,linePart3];
end

if any(strfind(line,'#rpA2#'))
    paramInd = strfind(line,'#rpA2#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.crackTipZoneSizes.rpA2,formatSpec);
    linePart3 = line((paramInd+6):end);
    line = [linePart1,linePart2,linePart3];
end

if any(strfind(line,'#rpB1#'))
    paramInd = strfind(line,'#rpB1#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.crackTipZoneSizes.rpB1,formatSpec);
    linePart3 = line((paramInd+6):end);
    line = [linePart1,linePart2,linePart3];
end

if any(strfind(line,'#rpB2#'))
    paramInd = strfind(line,'#rpB2#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.crackTipZoneSizes.rpB2,formatSpec);
    linePart3 = line((paramInd+6):end);
    line = [linePart1,linePart2,linePart3];
end

%Element sizes
if any(strfind(line,'#elSize1#'))
    error('#elSize1# found in master .py file. This has been superseded. Aborting.')
end

if any(strfind(line,'#elSize2#'))
    error('#elSize2# found in master .py file. This has been superseded. Aborting.')
end

if any(strfind(line,'#elSize3#'))
    paramInd = strfind(line,'#elSize3#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.elementSizes.elSize3,formatSpec);
    linePart3 = line((paramInd+9):end);
    line = [linePart1,linePart2,linePart3];
end

%Mesh transition positions
if any(strfind(line,'#meshTransA#'))
    error('#meshTransA# found in master .py file. This has been superseded by #meshTransSketchA#')
end

if any(strfind(line,'#meshTransB#'))
    error('#meshTransB# found in master .py file. This has been superseded by #meshTransSketchB#')
end

if any(strfind(line,'#meshTransSketchA#'))
    paramInd = strfind(line,'#meshTransSketchA#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.meshTransitionPositions.meshTransSketchA,formatSpec);
    linePart3 = line((paramInd+18):end);
    line = [linePart1,linePart2,linePart3];
end

if any(strfind(line,'#meshTransSketchB#'))
    paramInd = strfind(line,'#meshTransSketchB#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.meshTransitionPositions.meshTransSketchB,formatSpec);
    linePart3 = line((paramInd+18):end);
    line = [linePart1,linePart2,linePart3];
end

if any(strfind(line,'#meshTransZ#'))
    paramInd = strfind(line,'#meshTransZ#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.meshTransitionPositions.meshTransZ,formatSpec);
    linePart3 = line((paramInd+12):end);
    line = [linePart1,linePart2,linePart3];
end

%Expression for crack face pressure load
if any(strfind(line,'#CFPstring#'))
    error('#CFPstring# found in master .py file. This has been superseded. Aborting.')
end

%Poisson's ratio
if any(strfind(line,'#nu#'))
    paramInd = strfind(line,'#nu#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.materialParams.nu,formatSpec);
    linePart3 = line((paramInd+4):end);
    line = [linePart1,linePart2,linePart3];
end

%Plate sizes
if any(strfind(line,'#b#'))
    paramInd = strfind(line,'#b#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.plateSizes.b,formatSpec);
    linePart3 = line((paramInd+3):end);
    line = [linePart1,linePart2,linePart3];
end

if any(strfind(line,'#plEdgeA#'))
    paramInd = strfind(line,'#plEdgeA#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.plateSizes.plHalfWidth,formatSpec);
    linePart3 = line((paramInd+9):end);
    line = [linePart1,linePart2,linePart3];
end

if any(strfind(line,'#plEdgeB#'))
    paramInd = strfind(line,'#plEdgeB#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(-modelParamStruct.plateSizes.plHalfWidth,formatSpec);
    linePart3 = line((paramInd+9):end);
    line = [linePart1,linePart2,linePart3];
end

if any(strfind(line,'#plDepth#'))
    paramInd = strfind(line,'#plDepth#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.plateSizes.plDepth,formatSpec);
    linePart3 = line((paramInd+9):end);
    line = [linePart1,linePart2,linePart3];
end

if any(strfind(line,'#boundBoxMag#'))
    paramInd = strfind(line,'#boundBoxMag#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = [num2str(-modelParamStruct.plateSizes.boundBoxMag,formatSpec),',',...
        num2str(-modelParamStruct.plateSizes.boundBoxMag,formatSpec),',',...
        num2str(-modelParamStruct.plateSizes.boundBoxMag,formatSpec),',',...
        num2str(modelParamStruct.plateSizes.boundBoxMag,formatSpec),',',...
        num2str(modelParamStruct.plateSizes.boundBoxMag,formatSpec),',',...
        num2str(modelParamStruct.plateSizes.boundBoxMag,formatSpec)];
    linePart3 = line((paramInd+13):end);
    line = [linePart1,linePart2,linePart3];
end

%Number of contours in History Output Requests
if any(strfind(line,'#noContoursA#'))
    paramInd = strfind(line,'#noContoursA#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.outputRequests.noContoursA);
    linePart3 = line((paramInd+13):end);
    line = [linePart1,linePart2,linePart3];
end
if any(strfind(line,'#noContoursB#'))
    paramInd = strfind(line,'#noContoursB#');
    linePart1 = line(1:(paramInd-1));
    linePart2 = num2str(modelParamStruct.outputRequests.noContoursB);
    linePart3 = line((paramInd+13):end);
    line = [linePart1,linePart2,linePart3];
end

%Non-linear geometry
if any(strfind(line,'#nlgeomFlag#'))
    paramInd = strfind(line,'#nlgeomFlag#');
    linePart1 = line(1:(paramInd-1));
    if modelParamStruct.nlgeomFlag
        linePart2 = 'ON';
    else
        linePart2 = 'OFF';
    end
    linePart3 = line((paramInd+12):end);
    line = [linePart1,linePart2,linePart3];
end

%% .py file line for output
lineOut = line;

end