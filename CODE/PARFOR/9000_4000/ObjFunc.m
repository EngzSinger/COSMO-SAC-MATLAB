function [Err]=ObjFunc(C)


global combFlag parameterFlag cosmoFlag ExpDataPath Aeff;

if length(C)>2
	Aeff=C(3);
else
	Aeff=7.25;
end

combFlag=1;
parameterFlag=4;
cosmoFlag=3;
ExpDataPath='/home/zxwu/WORK_SPACE/COSMO_SAC/EXPDATA/IDAC/';
CationCosmoPath='/home/zxwu/WORK_SPACE/COSMO_SAC/DATABASE/COSMO/CATION/';
AnionCosmoPath='/home/zxwu/WORK_SPACE/COSMO_SAC/DATABASE/COSMO/ANION/';
SoluteCosmoPath='/home/zxwu/WORK_SPACE/COSMO_SAC/DATABASE/COSMO/SOLUTE/';
SystemIDAC=ReadIDACList('TrainList.txt');
%SystemIDAC=ReadIDACList('IDACList.txt');
Nsys=length(SystemIDAC{1,5});
Errg=zeros(Nsys,1);
disp(C);
if C(1)<100 | C(2)<100
	Err=1000;
else
Mypar=parpool('local',24);
parfor i=1:Nsys
	T=SystemIDAC{1,4}(i);
	%comp1=SystemIDAC{1,1}{i};
	%disp(comp1);
	%[area1,volume1,prf1]=PurePrf(comp1,SoluteCosmoPath,cosmoFlag);
	[area1,volume1,prf1]=PurePrf(SystemIDAC{1,1}{i},SoluteCosmoPath,cosmoFlag);
	%cation=SystemIDAC{1,2}{i};
	%anion=SystemIDAC{1,3}{i};
	%comp2=['[',cation,'][',anion,']'];
	%[Acation,Vcation,prf_cation]=PurePrf(cation,CationCosmoPath, cosmoFlag);
	[Acation,Vcation,prf_cation]=PurePrf(SystemIDAC{1,2}{i},CationCosmoPath, cosmoFlag);
	%[Aanion,Vanion,prf_anion]=PurePrf(anion,AnionCosmoPath,cosmoFlag);
	[Aanion,Vanion,prf_anion]=PurePrf(SystemIDAC{1,3}{i},AnionCosmoPath,cosmoFlag);
	area2=Acation+Aanion;
	volume2=Vcation+Vanion;
	prf2=prf_cation+prf_anion;
	x1=1e-6;
	X=[x1,1-x1];
	Acosmo=[area1(1),area2(1)];
	Vcosmo=[volume1,volume2];
	[gammaR,lngammaR]=gammares2007global(Vcosmo,Acosmo,X,prf1,prf2,T,parameterFlag,C);
	[gammaC,lngammaC]=gammacomb(Vcosmo,Acosmo,X,combFlag);
	lngamma=lngammaR+lngammaC;
    	gamma=exp(lngamma);
	Gexp=SystemIDAC{1,5}(i);
	g_Err=abs(gamma(1)-Gexp)/Gexp;
	Errg(i)=g_Err;
	%disp(cation);
	%disp(anion);
	%disp(comp1);
%	disp(Gexp);
%	disp(T);
%	disp(g_Err);
%	disp(gamma);
end

Err=sum(Errg)/Nsys;
delete(Mypar);
end
%parpool closed;
end
