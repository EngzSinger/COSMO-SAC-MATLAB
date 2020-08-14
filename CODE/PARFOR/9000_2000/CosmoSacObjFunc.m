function [Err]=CosmoSacObjFunc(C)
% �˺���������Ż���Ŀ�꺯��
% C�����洢COSMO-SACģ���໥���õ���������
% ������ȫŨ�ȷ�Χ��ȫ���ܵ���ϵ��Ϊ��ϵ����
% 
% clear;clc
% ȫ�ֱ���
global combFlag parameterFlag cosmoFlag ExpDataPath SoluteIndex Aeff;
%global combFlag parameterFlag cosmoFlag ExpDataPath SoluteIndex Aeff;

%+ Effective area of fragments
if length(C)>2
	Aeff=C(3);
else
	Aeff=7.25;
end

%+ѡ������������ϵ��Ĺ�ʽ
combFlag=1;
%+���㽻�����õĲ���ѡ��,parameterFlag(=1��COSMO-SAC(2007)�еĲ���=2��Wang Shu C++�����еĲ���=3�����粩ʿ���ĵĲ���)
parameterFlag=4;        %   =0�������Ż�
%+ѡ��cosmo�ļ���ԴcosmoFlag(=1��COSMOTherm�����=2��VT-2005��=3��G09���)
cosmoFlag=3;
%+ʵ������ļ���·��
%ExpDataPath='/home/whtu/COSMO-MATLAB/EXPDATA/GLE/';
ExpDataPath='/home/zxwu/WORK_SPACE/MATLAB_COSMO/EXPDATA/GLE/';
%+cosmo�ļ�·��
%CationCosmoPath='~/COSMO-MATLAB/COSMO/G09/Cation/6-311Gd/';
CationCosmoPath='/home/zxwu/WORK_SPACE/MATLAB_COSMO/DATABASE/COSMO/CATION/6-311Gd/';
%AnionCosmoPath='~/COSMO-MATLAB/COSMO/G09/Anion/6-311pGd/';
AnionCosmoPath='/home/zxwu/WORK_SPACE/MATLAB_COSMO/DATABASE/COSMO/ANION/6-311pGd/';
%SoluteCosmoPath='~/COSMO-MATLAB/COSMO/G09/Solute/6-311Gd/';
SoluteCosmoPath='/home/zxwu/WORK_SPACE/MATLAB_COSMO/DATABASE/COSMO/SOLUTE/6-311Gd/';
%+
SoluteIndex=ReadSoluteIndex('SoluteIndex.txt');     % SoluteIndex is a cell 1-by-SoluteNumber
% 
%   ����VLE
%����Ҫ�����ϵͳ
SystemVle=ReadSystemVleList('GLEList.txt');
Nsys=length(SystemVle{1,5});
ErrP=zeros(Nsys,1);
for i=1:Nsys  %����ÿ��ϵͳ���
    T=SystemVle{1,5}(i);  %ϵͳ�¶ȣ�K
	comp1=SystemVle{1,2}{i};	% Soute Name e.g. Methanol
	%+������1��profile
    [area1,volume1,prf1]=PurePrf(comp1,SoluteCosmoPath,cosmoFlag);
    % Psat=SystemVle{1,6}(i); %���ʱ�������ѹ��KPa
	[Psat,TLimit,Parameters]=antoine(comp1,T,SystemVle{1,6}(i));
    %+COmpound 2 is Ionic Liquid    
    cation=SystemVle{1,3}{i};     % Cation Abbreviation e.g. C4MIm
    anion=SystemVle{1,4}{i};      % Aation Abbreviation e.g. Tf2N
    comp2=['[',cation,'][',anion,']'];       %  Ionic Liquids  Abbreviation e.g. [C4MIm][Tf2N]
    %+�������Һ���profile
    [Acation,Vcation,prf_cation]=PurePrf(cation,CationCosmoPath, cosmoFlag);
    [Aanion,Vanion,prf_anion]=PurePrf(anion,AnionCosmoPath,cosmoFlag);
    area2=Acation+Aanion;
    volume2=Vcation+Vanion;
    prf2=prf_cation+prf_anion;
    %+����ʵ�����
    ExpPath=[ExpDataPath,SystemVle{1,1}{i},'.txt'];
    ExpData=load(ExpPath);
    % %+�������
    x1=ExpData(:,1);
    X=[x1,1-x1];
    x1Len=length(x1);
    DataNum=sum(x1>0);     %   ��ϵ��ݵ����,ȥ��x1=0.0�ĵ�
    deltLen=x1Len-DataNum;

    %����cosmo������
    Acosmo=[area1(1),area2(1)];
    Vcosmo=[volume1,volume2];
    %����gammares����
    % [gammaR,lngammaR]=gammares(Vcosmo,Acosmo,X,prf1,prf2,T);
    %����gammares2007��������gammares����ʵ����ͬ�Ĺ��ܣ�ֻ���ڱ�ʾ�໥����ʱ��ͬ
    [gammaR,lngammaR]=gammares2007global(Vcosmo,Acosmo,X,prf1,prf2,T,parameterFlag,C);
    %����gammacomb����
    [gammaC,lngammaC]=gammacomb(Vcosmo,Acosmo,X,combFlag); 
    lngamma=lngammaR+lngammaC;
    gamma=exp(lngamma);
    Pcal=x1.*gamma(:,1).*Psat;
    Pcal(find(Pcal>Psat))=Psat;
    p_Err=abs(Pcal((deltLen+1):x1Len)-ExpData((deltLen+1):x1Len,2))./ExpData((deltLen+1):x1Len,2).*100.0/DataNum;
    ErrP(i)=sum(p_Err);    
    
end

Err=sum(ErrP);
end
