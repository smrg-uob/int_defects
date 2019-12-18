function [ C ] = read_py( filename )
%read_py.m
%Harry Coules 2015
%
%DESCRIPTION
%This script is used to read the lines of a python script into Matlab. It
%returns the cell aray C which contains the lines of the .py file as
%strings.
%
%INPUT ARGUMENTS
% filename - Name of the .py file as a string
%
%OUTPUT ARGUMENTS
% C - Cell array of strings containing the file lines.

%% READ IN THE LINES OF THE .PY FILE
fid = fopen(filename,'r','n','UTF-8');
C = textscan(fid, '%s', Inf, 'Delimiter','\n'); C = C{1};
fclose(fid);

end