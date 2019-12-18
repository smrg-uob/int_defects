function int_defects_write_py( modelParamStruct, filenameMaster, filenameWrite)
%int_defects_write_py.m
%Harry Coules 2015
%
%DESCRIPTION
%This function is used to write an Abaqus python script using an existing
%'master' (i.e. template) script, and parameter defined in
%modelParamStruct. It writes the new file into the current directory.
%Certain strings eg. #aA1# in the master file indicate where the parameters
%should be written to.
%
%INPUT ARGUMENTS
% modelParamScript - Struct containing prameters to be written to the .py
%   file.
% filenameMaster - File name of the master python script (including file
%   extension).
% filenameWrite - File name of the python script to be written (including
%   the file extension).

%% Read in the master python script
[ C0 ] = read_py(filenameMaster);

%% Create new script from template line-by-line
C1 = {};
for k1 = 1:length(C0)
    line = char(C0(k1));
    
    %Make 'basic' changes to lines of the master .py file
    line = int_defects_py_line_basic(line, modelParamStruct);
    
    %Make more detailed changes to the lines of the master .py file. For
    %example, modifying the .py file to use getClosest rather than
    %getSequenceFromMask for selecting eges when creating mesh seeds.
    if isfield(modelParamStruct,'edgeCoords')
        line = int_defects_py_line_edges(line, modelParamStruct);
    end
    
    %Write modified line to new cell array
    C1 = [C1;line];
end

%% Write new python script
fid = fopen(filenameWrite,'w');
for k1 = 1:size(C1,1);
    fprintf(fid,[C1{k1},'\n']);
end
fclose(fid);

end