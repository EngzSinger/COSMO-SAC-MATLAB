function [C]=ReadSystemVleList(FileName)
% This function read the data as a format to a cell data C from FileName
fid=fopen(FileName);
if fid<0
    disp(['THE ReadSystemVleList FAILS TO OPEN THE FILE:  ',FileName]);
else
    C=textscan(fid,'%s %s %s %s %f %f');
end
fclose(fid);
end