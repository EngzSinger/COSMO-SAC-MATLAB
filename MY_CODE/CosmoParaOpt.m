%function be optimized with two parameter

% initial value for opt

ac0 = [9392.507,3631.760];

% read global parameter and distribute it to every function
% avoid repeatly open and close file
para = Para;

%file name of train list
ListFile = 'MeoList.txt';
fList = fopen(ListFile);
if fList < 0 
	error(['Please Check ',ListFile,' exist!!!']);
else
	trainList = textscan(fList,'%s %s %s %f %f %f');
end
fclose(fList);

% replace solute anion cation abbrevation to their index
% conviently find related cosmo file

%SoluteAbbr = unique(trainList{1,1});
%CationAbbr = unique(trainList{1,2});
%AnionAbbr = unique(trainList{1,3});

SolutePath = para('SolutePath');
CationPath = para('CationPath');
AnionPath = para('AnionPath');

% deal with solute abbr
fSolute = fopen([SolutePath,'Index.txt']);
if fSolute < 0
	error(['Please Check Index.txt exist in ',SolutePath]);
else
	SoluteCell = textscan(fSolute,'%s %s');
end
fclose(fSolute);

SoluteMap = containers.Map(SoluteCell{2},SoluteCell{1});
for i = 1:length(trainList{1,2})
	trainList{1,1}{i} = SoluteMap(trainList{1,1}{i});
end

% deal with cation abbr
fCation = fopen([CationPath,'Index.txt']);
if fCation < 0
        error(['Please Check Index.txt exist in ',CationPath]);
else
        CationCell = textscan(fCation,'%s %s');
end
fclose(fCation);

CationMap = containers.Map(CationCell{2},CationCell{1});
for i = 1:length(trainList{1,2})
        trainList{1,2}{i} = CationMap(trainList{1,2}{i});
end

%deal with anion abbr
fAnion = fopen([AnionPath,'Index.txt']);
if fAnion < 0
        error(['Please Check Index.txt exist in ',AnionPath]);
else
        AnionCell = textscan(fAnion,'%s %s');
end
fclose(fAnion);

AnionMap = containers.Map(AnionCell{2},AnionCell{1});
for i = 1:length(trainList{1,2})
        trainList{1,3}{i} = AnionMap(trainList{1,3}{i});
end


options = optimset;
options = optimset(options,'Display','iter','TolFun',1e-3,'TolX',1e-3);

tic;
disp('Optimizing ! PLEASE WAIT >>>>>>>>>>>>>>>>>');
[ac,fval,exitflag,output] = fminsearch(@(ac) OptFunc(ac,para,trainList),ac0,options)
toc
