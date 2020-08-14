function [gamma,lngamma]=gammares2007(Vcosmo,Acosmo,X,prf1,prf2,T,parameterFlag)
%��1������VcosmoΪ��ֵ����������Vcosmo(1)��Ӧ���1��cosmo���Vcosmo(2)��Ӧ���2��cosmo���
%��2������AcosmoΪ��ֵ����������Acosmo(1)��Ӧ���1��cosmo���Acosmo(2)��Ӧ���2��cosmo���
%��3������XΪ��ֵ���ɾ���ÿһ�д��һ����ɵ㣻
%��4������prf1��Ӧ�������1��sigma-profile��Ϊ2ά����column1��Ӧ����꣬column2Ϊ������profile��column3Ϊ���ֵ�profile����λΪ���λ��A^2
%��5������prf2��Ӧ�������2��sigma-profile��Ϊ2ά����column1��Ӧ����꣬column2Ϊ������profile��column3Ϊ���ֵ�profile����λΪ���λ��A^2
%��6������mixprf��Ӧ�������1�����2��ϵ�sigma-profile��Ϊ2ά����column1��Ӧ����꣬column2Ϊ������profile��column3Ϊ���ֵ�profile����λΪ���λ��A^2
%��7������TΪ��ϵ���¶ȣ���λΪK

global Aeff;
%constant for calculation
%Aeff=7.25;              %��Ч��׼Ƭ�������Ƭ��ƽ������ֵ����һ��,from the paper of COSMO-SAC(2007)
%% Information for compund and their profile
numComp=length(Vcosmo); %�����
numPoint=length(X(:,1));%��ݵ���
sumA=X*(Acosmo');       %��Ϻ������,����������Ӧÿ����Ļ�����
%% Rebuild the pure compund sigma-profile
[numSig,prfLevel]=size(prf1);
prfLevel=prfLevel-1;
numSig=prfLevel*numSig;
pureProfile=zeros(numSig,numComp);
pureProfile(:,1)=[prf1(:,2);prf1(:,3)];
pureProfile(:,2)=[prf2(:,2);prf2(:,3)];
%% Now caculating the residual activity coeffecient for each compund
% Initialize 
gamma=ones(numPoint,numComp);
lngamma=zeros(numPoint,numComp);
sumgamma=zeros(numPoint,numComp);
% Calcuating segmments activity coeffecient of pure compund when the
% terperature is constance for all data
if (length(T)==1)  
	[seggammaPure1]=seggamma(prf1,Acosmo(1),T,parameterFlag);
	[seggammaPure2]=seggamma(prf2,Acosmo(2),T,parameterFlag);    
    seggammaPure=[seggammaPure1;seggammaPure2];    
end
% xlswrite('seggammaPure_new.xls',seggammaPure);
for i=1:numPoint
    % Calcuating segmments activity coeffecient of pure compund
    if (length(T)>1)
		[seggammaPure1]=seggamma(prf1,Acosmo(1),T(i),parameterFlag);
		[seggammaPure2]=seggamma(prf2,Acosmo(2),T(i),parameterFlag);
		seggammaPure=[seggammaPure1;seggammaPure2];                
    end
    % Calcuating segmments activity coeffecient of mixture
    mixprf=prf1*X(i,1)+prf2*X(i,2);   % mixture sigma-profile
    if(length(T)>1)
        [seggammaMix]=seggamma(mixprf,sumA(i),T(i),parameterFlag);
    else
        [seggammaMix]=seggamma(mixprf,sumA(i),T,parameterFlag);
    end
    % Calculating residual activity coeffecient of each compund
    for k=1:numComp
        for j=1:numSig
            sumgamma(i,k)=sumgamma(i,k)+((pureProfile(j,k)/Aeff)*(log(seggammaMix(j)/seggammaPure(k,j))));
        end
        gamma(i,k)=exp(sumgamma(i,k));
        lngamma(i,k)=log(gamma(i,k));
    end
    
end

end
    %���ϵ��ʣ����������
    
        
       