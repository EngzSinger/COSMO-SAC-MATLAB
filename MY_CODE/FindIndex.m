function [Index] = FindIndex(Path,CompName)
% Return the comp's index in the path
% Index.file should contain following information
% compname	Index
% PF6    	anion0001

searchFlag = 0;
fid = fopen([Path,'Index.txt']);
if fid<0
	error(['Please check whether exist Index.txt in ',Path]);
end
while ~feof(fid)
	C = textscan(fid,'%s %s');
end
for i = 1:length(C{1,1})
	if strcmp(C{1,2}{i},CompName)
		Index = C{1,1}{i};
		searchFlag = 1;
		break;
	end
end
if searchFlag<1
	error(['Please be sure ' ,CompName,'in the file']);
end
fclose(fid);
end
