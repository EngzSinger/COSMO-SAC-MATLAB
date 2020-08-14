function [atom,segments,Vcosmo,Acosmo,atomNum,segNum] = ReadG09(CompName,Path,para)
% This function read Comp's atom information,segments for further calculation

%para = Para();
u2A = str2num(para('u2A'));

AtomList = {'H','B','C','N','O','F','Al','P','S','Cl','Fe','Zn','Ga','In'};
AtomIndex = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};
AtomMap = containers.Map(AtomList,AtomIndex);
%[Index] = FindIndex(Path,CompName);
CosmoPath = [Path,CompName,'.cosmo'];
fid = fopen(CosmoPath);

%check if file of comp exist
if(fid<0)
	error(['Please Check whether exist',CompName,'.cosmo(',CompName ,'.cosmo) File in' Path]);
end

% read account of atom and segments
while ~feof(fid)
	tline = fgetl(fid);
	if strfind(tline,'!DATE')
		AtomCount = 0;
		while ~feof(fid);
			tline = fgetl(fid);
			if strfind(tline,'end')
				break;
			end
			AtomCount = AtomCount+1;
		end
	end

	if strfind(tline,'length')
		for i = 1:4
			tline = fgetl(fid);
		end
		SegCount = 0;
		while ~feof(fid)
			tline = fgetl(fid);
			SegCount = SegCount+1;
		end
		break
	end
end

atomNum = AtomCount;
segNum = SegCount;


% read cosmo file again to get atom information and segments
fseek(fid,0,-1);
atom = zeros(atomNum,4);
segments = zeros(segNum,9);
while ~feof(fid)
	tline = fgetl(fid);
	if strfind(tline,'area  =')
		[str1,str2,Acosmo] = strread(tline,'%s%s%f');
		Acosmo = u2A * u2A * Acosmo;
	end
	if strfind(tline,'volume=')
		[str1,Vcosmo] = strread(tline,'%s%f');
		Vcosmo = u2A * u2A * u2A * Vcosmo;
	end
	if strfind(tline,'!DATE')
		[C] = textscan(fid,'%s %f %f %f %s %d %s %s %f',atomNum);
		atom(:,2) = C{1,2};
		atom(:,3) = C{1,3};
		atom(:,4) = C{1,4};
		str = C{1,8};
		for j = 1:atomNum
			str4 = char(str(j));
			atom(j,1) = AtomMap(str4);
		end
	end
	if strfind(tline,'length')
		for i = 1:4
			tline = fgetl(fid);
		end
		[C] = textscan(fid,'%d %d %f %f %f %f %f %f %f',segNum);
		for i = 1:9
			segments(:,i) = C{1,i};
		end
		break
	end
end
fclose(fid);
end
