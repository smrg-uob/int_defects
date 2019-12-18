function [ line2 ] = int_defects_py_line_builder( coordArray, varargin )
%int_defects_py_line_builder.m
%Harry Coules 2015
%
%DESCRIPTION
%This function is used to create a string containing coordinates, deformation
%plasticity parameters or plastic flow data in a certain format used by the
%Abaqus scripting interface
%
%INPUT ARGUMENTS
%   coordArray - Array of coordinates to be written to the .py file in
%       standard format.
%   OPTIONAL INPUT ARGUMENTS
%   formatSpec - String specifying the number format which should be used.
%       If this omitted, the default is '%.7f'.
%   tableTypeStr - String indicating what sort of data table is being
%       processed:
%       'coords'   =  n-by-3 table of spatial coordinates (3D).
%       'defplast' =  Column vector of parameters defining a deformation plasticity model.
%       'incplast' =  n-by-2 table of flow stress vs. strain data for an incremental plasticity model.
%
%OUTPUT ARGUMENTS
%   line2 - String containing the line of the .py file generated from
%       coordArray.
%
%% Optional input arguments
if ~isempty(varargin)
    if length(varargin) == 1
        formatSpec = varargin{1};
        tableTypeStr = 'coords';
    elseif length(varargin) == 2
        formatSpec = varargin{1};
        tableTypeStr = varargin{2};
    else
        error('Too many input arguments.')
    end
else
    formatSpec = '%.7f';    %Silently default to this if the optional input argument is omitted.
    tableTypeStr = 'coords';
end

%% Generate the .py file line
if strcmpi(tableTypeStr,'incplast')
    line2 = '(';
    for k1 = 1:size(coordArray,1)
        line2 = [line2,'('];
        for k2 = 1:size(coordArray,2)
            line2 = [line2,num2str(coordArray(k1,k2),formatSpec)];
            if k2==1
                line2 = [line2,','];
            end
        end
        line2 = [line2,'),'];
    end
    line2 = [line2,'))'];
elseif strcmpi(tableTypeStr,'defplast')
    line2 = '((';
    for k1 = 1:size(coordArray,1)
        line2 = [line2,num2str(coordArray(k1,1),formatSpec)];
        if k1~=size(coordArray,1)
            line2 = [line2,','];
        end
    end
    line2 = [line2,'),))'];
elseif strcmpi(tableTypeStr,'coords')
    line2 = '';
    for k1 = 1:size(coordArray,1)
        line2 = [line2,'(('];
        for k2 = 1:size(coordArray,2)
            line2 = [line2,num2str(coordArray(k1,k2),formatSpec)];
            if k2~=3
                line2 = [line2,','];
            end
        end
        line2 = [line2,'),)'];
        if k1~=size(coordArray,1)
            line2 = [line2,','];
        end
    end
    line2 = [line2,')'];
else
    error('Unrecognised input table type specified.')
end

end