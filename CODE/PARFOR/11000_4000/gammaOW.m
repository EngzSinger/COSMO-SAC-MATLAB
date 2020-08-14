function [gamma]=gammaOW(comp3,x)
% This function aims to calculate the Octanol-Water partition coefficient
% of comp3 at 298.15 K
% x(1,1),x(1,2),x(1,3) are the mole fraction of octano,water,comp3 in octanol-reach phase,respectively
% x(2,1),x(2,2),x(2,3) are the mole fraction of octano,water,comp3 in
% octanol-reach phase,respectively
%+ѡ������������ϵ���Ĺ�ʽ combFlag=1;
%+��ϵ����ѡ��ILflag��=1������Һ����ϵ��=0��������Һ����ϵ��ILflag=1;
%+���㽻�����õĲ���ѡ��,parameterFlag(=1��COSMO-SAC(2007)�еĲ�����=2��Wang Shu C++�����еĲ�����=3�����粩ʿ���ĵĲ���) parameterFlag=2;
%+ѡ��cosmo�ļ���ԴcosmoFlag(=1��COSMOTherm�����=2��VT-2005��=3��G09����) cosmoFlag=2;
global combFlag parameterFlag cosmoFlag Aeff;
%constant for calculation
%Aeff=7.25;              %��Ч��׼Ƭ���������Ƭ��ƽ��������ֵ����һ��,from the paper of COSMO-SAC(2007)
rgas=0.001987;          %���峣����kcal/mol/K
tol=0.000000001;            %������׼
numComp=3;      %�����
numPoint=2;     %���ݵ�
T=298.15;
comp1='octanol';
comp2='water';
[area1,volume1,prf1]=PurePrf(comp1,cosmoFlag);
[area2,volume2,prf2]=PurePrf(comp2,cosmoFlag);
[area3,volume3,prf3]=PurePrf(comp3,cosmoFlag);
Acosmo=[area1(1),area2(1),area3(1)];
Vcosmo=[volume1,volume2,volume3];
sumA=x*(Acosmo');       %��Ϻ������,����������Ӧÿ����Ļ�����
sigma=prf1(:,1);        %sigma-profile����������
prfLevel=2;             %profile��ֳɵ�ά�ȣ�����ȡ2��ʾprofile��ֳɷ������������֣�N O F�Լ�����������H��
numSig=prfLevel*length(sigma);   %��profile���������ݵ���
%
%
%���㴿���Ƭ�εĻ��ϵ��
%     if i==1
    convPure=ones(numComp,numSig);
    seggammaPureOld=ones(numComp,numSig);
    %+���촿��ֵ�sigma-profile����
    purePrf=zeros(numSig,numComp);  %��λ��ȻΪ�����λ��A^2
    purePrf(:,1)=[prf1(:,2);prf1(:,3)];
    purePrf(:,2)=[prf2(:,2);prf2(:,3)];
    purePrf(:,3)=[prf3(:,2);prf3(:,3)];
    seggammaPure=ones(numComp,numSig);
    for j=1:numComp        
%         seggammaPure=ones(numComp,numSig);
        istep=0;  %��¼��������
        while (sum(convPure(j,:))>=tol)
            istep=istep+1;
            convPure(j,:)=zeros(numSig,1);    %�б�Ҫ����Ϊ0��
            seggammaPureOld(j,:)=seggammaPure(j,:);
            for m=1:numSig
                summation=0.0;
                    for n=1:numSig
                        if (purePrf(n,j)~=0)
                            summation=summation+(purePrf(n,j)/Acosmo(j))*seggammaPureOld(j,n)*exp(-deltaW(m,n,sigma,parameterFlag)/rgas/T);  %Ҫ��鵥λ
                        end                            
                    end
                    seggammaPure(j,m)=exp(-log(summation));
                    seggammaPure(j,m)=(seggammaPure(j,m)+seggammaPureOld(j,m))/2.0;
            end
             convPure(j,:)=(seggammaPure(j,:)-seggammaPureOld(j,:)).^2;  %by Tu 2016-04-28
         end
            %convPure(j,:)=abs(seggammaPure(j,:)-seggammaPureOld(j,:))./seggammaPureOld(j,:);
           %convPure(j,:)=(seggammaPure(j,:)-seggammaPureOld(j,:)).^2;  %by Tu 2016-04-28
     end
    
%     end
    %���㴿���Ƭ�λ��ϵ������
%
gammaR=ones(numPoint,numComp);
lngammaR=zeros(numPoint,numComp);
sumgamma=zeros(numPoint,numComp);
for i=1:numPoint
    %������profile
    %+������profile��С
     TempMixprf=zeros(numSig,1);
    % mix_area=area1*X(1)+area2*X(2);%����������[�������������������������������]
    mixprf=prf1(:,2:3)*x(i,1)+prf2(:,2:3)*x(i,2)+prf3(:,2:3)*x(i,3);
    %������ʱ���prf
    TempMixprf=[mixprf(:,1);mixprf(:,2)];
    %������Ƭ�λ��ϵ��
    seggamma=ones(1,numSig);
    converge=ones(1,numSig);
    while sum(converge)>=tol
        seggammaold=seggamma;        
            for m=1:numSig
                summation=0.0;   %Ҫ���Ƿ������Ƿ����                
                    for n=1:numSig
                        if (mixprf(n)~=0)
                            summation=summation+mixprf(n)/sumA(i)*seggammaold(n)*exp(-deltaW(m,n,sigma,parameterFlag)/rgas/T);
                        end
                    end                
                seggamma(m)=exp(-log(summation));
                seggamma(m)=(seggamma(m)+seggammaold(m))/2.0;
            end        
        converge=(seggamma-seggammaold).^2;
    end
    %������Ƭ�λ��ϵ������
%
    %������ϵ��ʣ����
%     gamma=ones(1,numComp);
%     lngamma=ones(1,numComp);
%     sumgamma=ones(1,numComp);
    for k=1:numComp        
            for j=1:numSig
                sumgamma(i,k)=sumgamma(i,k)+((purePrf(j,k)/Aeff)*(log(seggamma(j)/seggammaPure(k,j)))); %ע���������ʽ�Ƿ���ȷ
            end        
       gammaR(i,k)=exp(sumgamma(i,k));  %         gamma(i,k)=exp(sumgamma(i,k)+lngamma(i,k));��������ʽ����ȷ��
       lngammaR(i,k)=log(gammaR(i,k));
    end
end
    %���ϵ��ʣ����������
%     
% ++��ʼ�����������ϵ��
gammaC=ones(numPoint,numComp);
lngammaC=zeros(numPoint,numComp);
for i=1:numPoint
[gammaC(i,:),lngammaC(i,:)]=gammacomb(Vcosmo,Acosmo,x(i,:),combFlag);
end

% ������ϵ��
lngamma=lngammaR+lngammaC;
gamma=exp(lngamma);

end