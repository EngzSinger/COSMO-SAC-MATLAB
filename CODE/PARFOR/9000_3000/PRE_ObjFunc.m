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
T=zeros(Nsys,1);
area1=zeros(Nsys,3);
volume1=zeros(Nsys,1);
prf1=zeros(Nsys,71,3);
Acation=zeros(Nsys,3);
Vcation=zeros(Nsys,1);
prf_cation=zeros(Nsys,71,3);
Aanion=zeros(Nsys,3);
Vanion=zeros(Nsys,1);
prf_anion=zeros(Nsys,71,3);
area2=zeros(Nsys,3);
volume2=zeros(Nsys,1);
prf2=zeros(Nsys,71,3);
x1=zeros(Nsys,1);
X=zeros(Nsys,2);
Acosmo=zeros(Nsys,2);
Vcosmo=zeros(Nsys,2);
gammaR=zeros(Nsys,2);
lngammaR=zeros(Nsys,2);
gammaC=zeros(Nsys,2);
lngammaC=zeros(Nsys,2);
gamma=zeros(Nsys,2);
Gexp=zeros(Nsys,1);

parfor i=1:Nsys
	T(i)=SystemIDAC{1,4}(i);
%	comp1=SystemIDAC{1,1}{i};
	[area1(i),volume1(i),prf1(i)]=PurePrf(SystemIDAC{1,1}{i},SoluteCosmoPath,cosmoFlag);
%	[area1,volume1,prf1]=PurePrf(comp1,SoluteCosmoPath,cosmoFlag);

%	cation=SystemIDAC{1,2}{i};
%	anion=SystemIDAC{1,3}{i};
%	comp2=['[',cation,'][',anion,']'];
%	[Acation,Vcation,prf_cation]=PurePrf(cation,CationCosmoPath, cosmoFlag);
	[Acation(i),Vcation(i),prf_cation(i)]=PurePrf(SystemIDAC{1,2}{i},CationCosmoPath, cosmoFlag);
%	[Aanion,Vanion,prf_anion]=PurePrf(anion,AnionCosmoPath,cosmoFlag);
	[Aanion(i),Vanion(i),prf_anion(i)]=PurePrf(SystemIDAC{1,3}{i},AnionCosmoPath,cosmoFlag);
	area2(i)=Acation(i)+Aanion(i);
	volume2(i)=Vcation(i)+Vanion(i);
	prf2(i)=prf_cation(i)+prf_anion(i);
	x1(i)=1e-6;
	X(i)=[x1,1-x1];
	Acosmo(i)=[area1(i,1),area2(i,1)];
	Vcosmo(i)=[volume1(i),volume2(i)];
	[gammaR(i),lngammaR(i)]=gammares2007global(Vcosmo(i),Acosmo(i),X(i),prf1(i),prf2(i),T(i),parameterFlag,C);
	[gammaC(i),lngammaC(i)]=gammacomb(Vcosmo(i),Acosmo(i),X(i),combFlag);
	lngamma(i)=lngammaR(i)+lngammaC(i);
    	gamma(i)=exp(lngamma(i));
	Gexp(i)=SystemIDAC{1,5}(i);
	Errg(i)=abs(gamma(i,1)-Gexp(i))/Gexp(i);
	%disp(cation);
	%disp(anion);
	%disp(comp1);
%	disp(Gexp);
%	disp(T);
%	disp(g_Err);
%	disp(gamma);
end

Err=sum(Errg)/Nsys;
end
