function exitStatus = int_defects_flip_normal_in_inp( fileToChange, normalsToFlip, contourName1, varargin)
%int_defects_flip_normal_in_inp.m
%Harry Coules 2015
%
%DESCRIPTION
%This function is used to flip the direction of the normal vector in the
%.inp file for contour integral evaluation at a crack or cracks.
%
%INPUT ARGUMENTS
%   fileToChange - Name of the .inp file to be modified (including
%       filename).
%   normalsToFlip - Logical of length 2 to define which cracks are to have
%       their normals flipped.
%   contourName1 - String defining the name of the first contour integral
%       output request.
%   *OPTIONAL*
%   contourName2 - String defining the name of the second contour integral
%       output request (for a two-crack model).
%
%OUTPUT ARGUMENTS
%   exitStatus - Logical to indicate whether the operation was successful
%   (true) or unsuccessful (false).
%
%NOTES
% - This function will only work for .inp files generated using the
%   int_defects software package. It relies upon there being exactly one or
%   two cracks in the model, and these having specific names.
% - If normalsToFlip is given as [1,0] or [0,1] then contourName1 is the
%   name of the appropriate crack's contour output request. If
%   normalsToFlip is given as [1,1] then contourName1 and contourName2
%   should be given.
%
%% Sanitise input arguments
%normalsToFlip
if length(normalsToFlip) == 1
    warning('normalsToFlip has a length of 1. Assuming that this refers to Crack A.')
    normalsToFlip = [normalsToFlip,false];
elseif length(normalsToFlip) ~= 2
    error('The input argument normalsToFlip should have a length of 1 or 2');
end

%Optional input argument: contourName2
if isempty(varargin)
    contourName2 = '';
elseif length(varargin)==1
    contourName2 = varargin{1};
else
    error('Too many input arguments provided.');
end

%% Contour strings
if normalsToFlip(1) && normalsToFlip(2)
    contourStrA = strcat('*Contour Integral, crack name=',contourName1,'_CrackA');
    contourStrB = strcat('*Contour Integral, crack name=',contourName2,'_CrackB');
elseif normalsToFlip(1)
    contourStrA = strcat('*Contour Integral, crack name=',contourName1,'_CrackA');
    contourStrB = '';
elseif normalsToFlip(2)
    contourStrA = '';
    contourStrB = strcat('*Contour Integral, crack name=',contourName2,'_CrackB');
else
    warning('normalsToFlip = [0,0]. Not attempting to modify the input file.');
    contourStrA = '';
    contourStrB = '';
end
exitStatus = true;

%% Read file
fid = fopen(fileToChange,'r','n','UTF-8');
C = textscan(fid, '%s', Inf, 'Delimiter','\n'); C = C{1};
fclose(fid);

%% Find and replace the crack normals in the *Contour integral keyword
cracksFlippedA = 0;
cracksFlippedB = 0;
for k1 = 1:length(C)
    %Find and replace normal for Crack A (if required)
    if normalsToFlip(1)
        if ~isempty(strfind(char(C(k1)),contourStrA))
            if ~isempty(strfind(char(C(k1+1)),'0., 0., 1.'))
                C{k1+1} = '0., 0., -1.';
                cracksFlippedA = cracksFlippedA + 1;
            elseif ~isempty(strfind(char(C(k1+1)),'0., 0., -1.'))
                C{k1+1} = '0., 0., 1.';
                cracksFlippedA = cracksFlippedA + 1;
            else
                warning('A *Contour Integral keyword for Crack A was found, but the normal was not given or an unexpected vector.')                
                exitStatus = false;
            end
        end
    end
    %Find and replace normal for Crack B (if required)
    if normalsToFlip(2)
        if ~isempty(strfind(char(C(k1)),contourStrB))
            if ~isempty(strfind(char(C(k1+1)),'0., 0., 1.'))
                C{k1+1} = '0., 0., -1.';
                cracksFlippedB = cracksFlippedB + 1;
            elseif ~isempty(strfind(char(C(k1+1)),'0., 0., -1.'))
                C{k1+1} = '0., 0., 1.';
                cracksFlippedB = cracksFlippedB + 1;
            else
                warning('A *Contour Integral keyword for Crack B was found, but the normal was not given or an unexpected vector.')   
                exitStatus = false;
            end
        end
    end
end

%% Check that if there has been an instruction to flip the normal for a crack, then at least one *Contour Integral for it has been changed
if normalsToFlip(1)
    if cracksFlippedA == 0
        warning('Not able to flip the normal vector for Crack A')
        exitStatus = false;
    end
end

if normalsToFlip(2)
    if cracksFlippedB == 0
        warning('Not able to flip the normal vector for Crack B')
        exitStatus = false;
    end
end

%% Write file
fid = fopen(fileToChange,'w');
for k1 = 1:size(C,1);
    fprintf(fid,[C{k1},'\n']);
end
fclose(fid);

end