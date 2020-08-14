function List=DataWash(ListFile)
% Data collected from variety resource can be significately deviate from each other
% A scientific and effecient method to wash experimental data is thus of great important
% In this method, the very data point will be treated in the consideration of their deviation
% while solute, anion, cation must agree to be the same data point, and of temperature can not greater than 5K

%ListFile='PF6List.txt';
lT=0.004;
hT=0.01;
lG=0.5;
hG=0.8;
LogPath='/home/zxwu/WORK_SPACE/COSMO_SAC/LOGFILE/';
LogFile=[ListFile,'.log'];
Log=[LogPath,LogFile];
fid=fopen(ListFile);
if fid<0
	disp(['FAILS TO OPEN THE FILE: ',ListFile]);
else
	List=textscan(fid,'%s %s %s %f %f %d %d');
end
fclose(fid);

% set Flag for every data point
NData=length(List{5});
DataFlag=0;
for i=1:NData
	Solute=List{1}{i};
	Cation=List{2}{i};
	Anion=List{3}{i};
	T=List{4}(i);
	gamma=List{5}(i);
	Flag=List{6}(i);
	if Flag ~= 0
		continue
	end
	DataFlag=DataFlag+1;
	List{6}(i)=DataFlag;
	for j=i+1:NData
		Solutej=List{1}{j};
		Cationj=List{2}{j};
        	Anionj=List{3}{j};
        	Tj=List{4}(j);
        	gammaj=List{5}(j);
        	Flagj=List{6}(j);
		if Flagj ~= 0
			continue
		end
		if strcmp(Solute,Solutej) & strcmp(Cation,Cationj) & strcmp(Anion,Anionj) & abs(T-Tj)<5
		List{6}(j)=DataFlag;
		end
	end
end

% check same point
for i=1:DataFlag
	Points=[find(List{6}(:)==i)];
	NPoint=length(Points);
	if NPoint == 1
		List{7}(Points)=1;
		continue
	end
	T=List{4}(Points);
	Td=AADEV(T);
	Gamma=List{5}(Points);
	Gd=AADEV(Gamma);
	if (Td < lT & Gd < lG) | (Td < hT & Td > lT & Gd < hG)
		List{4}(Points(1))=mean(T);
		List{5}(Points(1))=mean(Gamma);
		List{7}(Points(1))=1;
	end
end

RemovedIndex=find(List{7}(:)==0);

fLog=fopen(Log,'w');
fprintf(fLog,'%s\n','Abandoned Data Points Are As Follow: ');
fprintf(fLog,'%s\t%s\t%s\t%s\t%s\t%s\t\n','Solute','Cation','Anion','T','Gamma','DataFlag');
for i=1:length(RemovedIndex)
	fprintf(fLog,'%s\t%s\t%s\t%.2f\t%.2f\t%d\t\n',List{1}{RemovedIndex(i)},List{2}{RemovedIndex(i)},List{3}{RemovedIndex(i)},List{4}(RemovedIndex(i)),List{5}(RemovedIndex(i)),List{6}(RemovedIndex(i)));
end
fclose(fLog);
end
% read result data
%RIndex=find(List{7}(:)==1);
%RList=List(RIndex);

%function ad = AADEV(a)

%	ave = mean(a)*ones(length(a),1)
%	dev = abs(a-ave)
%	ad = sum(dev)/mean(a)
%end
