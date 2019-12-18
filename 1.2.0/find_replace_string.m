function [varargout] = find_replace_string( oldStr, newStr, fileToChange )
%find_replace_string.m
%Harry Coules 2015
%
%DESCRIPTION
%   This function performs a find-and-replace operation within a named text
%   file.
%
%INPUT ARGUMENTS
%   oldStr - The string to be replaced.
%   newStr - The new string that the old one will be replaced with.
%   fileToChange - Filename of the file to be modified (including the file
%   extension).
%
%OUTPUT ARGUMENTS
%   OPTIONAL
%   noReplacements - The number of lines in the file in which oldStr was found.

%% Read file
fid = fopen(fileToChange,'r','n','UTF-8');
C0 = textscan(fid, '%s', Inf, 'Delimiter','\n'); C0 = C0{1};
fclose(fid);

%% Find and replace string
C1 = C0;
replacedFlag = zeros(length(C0),1);
for k1 = 1:length(C0)
    if  any(strfind(char(C0(k1)), oldStr))
        replacedFlag(k1) = 1;
        C1{k1} = strrep(char(C0(k1)), oldStr, newStr);
    end
end

%% Modify file
fid = fopen(fileToChange,'w');
for k1 = 1:length(C1)
    fprintf(fid,[C1{k1},'\n']);
end
fclose(fid);

%% If requested, provide the number of instances of the search string
if nargout == 1
    varargout{1} = sum(replacedFlag);
end

end

