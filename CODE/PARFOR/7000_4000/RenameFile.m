function []=RenameFile(RenameList,FileType)
% This function renames the files named based on first column to thenames of second column in RenameList file
% The format of the file RenameList as follows:
%                                               anion001	Cl
%                                               anion006	BF4
% That is to say the file anion001.txt is renamed Cl.txt
fid=fopen(RenameList);
NameIndex=textscan(fid,'%s %s');
fclose(fid);
for i=1:length(NameIndex{1,1})
    copyfile([NameIndex{1,1}{i},FileType],[NameIndex{1,2}{i},FileType]);
end
