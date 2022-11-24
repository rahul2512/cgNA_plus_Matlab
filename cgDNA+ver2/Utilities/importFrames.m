function [base, pho] = importFrames(filename)


formatSpec ='%*6s%12f%12f%12f%12f%12f%12f%12f%12f%12f%12f%12f%f%[^\n\r]';
%formatSpec = '%2f%4f%12f%12f%12f%12f%12f%12f%12f%12f%12f%12f%12f%f%[^\n\r]';

%% Open the text file.
fileID = fopen([ filename '.fra' ],'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

base = [dataArray{1:end-1}];

%% Open the text file.
fileID = fopen([ filename '.pfra' ],'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

pho = [dataArray{1:end-1}];

end

