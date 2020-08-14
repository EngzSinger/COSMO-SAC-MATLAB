%function [Err]=ObjFunc(C)


global combFlag parameterFlag cosmoFlag ExpDataPath Aeff;

%if length(C)>2
%	Aeff=C(3);
%else
%	Aeff=7.25;
%end
Aeff=7.25;
C=1;
combFlag=1;
parameterFlag=3;
cosmoFlag=3;
%ExpDataPath='/home/zxwu/WORK_SPACE/COSMO_SAC/EXPDATA/IDAC/';
CationCosmoPath='/home/zxwu/WORK_SPACE/COSMO_SAC/DATABASE/COSMO/CATION/';
AnionCosmoPath='/home/zxwu/WORK_SPACE/COSMO_SAC/DATABASE/COSMO/ANION/';
SoluteCosmoPath='/home/zxwu/WORK_SPACE/COSMO_SAC/DATABASE/COSMO/SOLUTE/';
SystemIDAC=ReadIDACList('TrainList.txt');
%SystemIDAC=ReadIDACList('tmp');
Nsys=length(SystemIDAC{1,5});
Errg=zeros(Nsys,1);
gam=zeros(Nsys,1);
disp(C);
for i=1:1
	T=SystemIDAC{1,4}(i);
	comp1=SystemIDAC{1,1}{i};
	[area1,volume1,prf1]=PurePrf(comp1,SoluteCosmoPath,cosmoFlag);

	cation=SystemIDAC{1,2}{i};
	anion=SystemIDAC{1,3}{i};
	comp2=['[',cation,'][',anion,']'];
	[Acation,Vcation,prf_cation]=PurePrf(cation,CationCosmoPath, cosmoFlag);
	[Aanion,Vanion,prf_anion]=PurePrf(anion,AnionCosmoPath,cosmoFlag);
	area2=Acation+Aanion;
	volume2=Vcation+Vanion;
	prf2=prf_cation+prf_anion;
	x1=1e-6;
	X=[x1,1-x1];
	Acosmo=[area1(1),area2(1)];
	Vcosmo=[volume1,volume2];
	[gammaR,lngammaR]=gammares2007global(Vcosmo,Acosmo,X,prf1,prf2,T,parameterFlag,C);
	[gammaC,lngammaC]=gammacomb(Vcosmo,Acosmo,X,0);
	lngamma=lngammaR+lngammaC;
    	gamma=exp(lngamma);
	Gexp=SystemIDAC{1,5}(i);
	gam(i)=gamma(1);
	g_Err=abs(gamma(1)-Gexp)/Gexp;
	Errg(i)=g_Err;
	%disp(cation);
	%disp(anion);
	%disp(comp1);
	disp(Gexp);
	disp(T);
	disp(g_Err);
	disp(gamma);
end

Err=sum(Errg);
	pathResult=['/home/zxwu/WORK_SPACE/COSMO_SAC/RESULT/Result_1.txt'];
	fid=fopen(pathResult,'w');
	fprintf(fid,'%s\n\n','Results summary:');
	fprintf(fid,'%s%d%s%d%s%s\n\n','Program contral parameterFlag=',parameterFlag,' ; cosmoFlag=',cosmoFlag,'     Calculated on: ',date);
	fprintf(fid,'%s\t%.4f\n\n','AARD=',Err/Nsys);
	fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n','No','Cation','Anion','Solute','T','Vcal','Vexp','ARD');
	for i=1:Nsys
		fprintf(fid,'%d\t%s\t%s\t%s\t%8.2f\t%.4f\t%.4f\t%.4f\n',i,SystemIDAC{1,2}{i},SystemIDAC{1,3}{i},SystemIDAC{1,1}{i},SystemIDAC{1,4}(i),gam(i),SystemIDAC{1,5}(i),Errg(i));
	end
	fclose(fid);
