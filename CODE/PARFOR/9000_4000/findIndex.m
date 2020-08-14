function [Index]=findIndex(Path,CompName)
% RETURN the CompName's index in Path
% e.g CompName          Index
%       C2MPip               cation284
%       TfO                     anion039
searchFlag=0;
fid=fopen([Path,'Index.txt']);
if fid<0
    error(['PLEASE CHECK WHETHER EXIST Index.txt FILE IN ',Path]);
end
while ~feof(fid)
    C=textscan(fid,'%s %s');
end
for i=1:length(C{1,1})
    if strcmp(C{1,2}{i},CompName)
        Index=C{1,1}{i};
        searchFlag=i;
        break;
    end
end
if searchFlag<1
    error(['CAN NOT FIND THE ',CompName,' IN ',Path,'Index.txt !']);
end
fclose(fid);
end