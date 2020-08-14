function [SoluteIndex]=ReadSoluteIndex(FileName)
% This function read the Index information and Antoine Parameters of Solute
% from FileName txt file
% The format of FileName txt file as follow£º
%       Total-Number-of-Solutes-is 36
%       Index               Name                    LowTemp     UpTemp      A(bar)      B(bar.K)        C(K)
%       Solute0001      Cyclohexane           293.06       354.73      3.96988     1203.526        -50.287
%       Solute0002      Benzene                  287.7        354.07      4.01814     1203.835        -53.226


fid=fopen(FileName);                               % Open the file
%while ~feof(fid)
    SoluteNum=textscan(fid,'%s %d',1);      % Get the Number of Solutes
%    tline=fgetl(fid);
    SoluteIndex=textscan(fid,'%s %s %f %f %f %f %f',SoluteNum{1,2}+1);    % Read the Information
%end 
fclose(fid);
