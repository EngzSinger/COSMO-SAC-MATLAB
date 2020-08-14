function [gamma,lngamma]=gammares2007global(Vcosmo,Acosmo,X,prf1,prf2,T,parameterFlag,C)
%��1������VcosmoΪ��ֵ����������Vcosmo(1)��Ӧ���1��cosmo���Vcosmo(2)��Ӧ���2��cosmo���
%��2������AcosmoΪ��ֵ����������Acosmo(1)��Ӧ���1��cosmo���Acosmo(2)��Ӧ���2��cosmo���
%��3������XΪ��ֵ���ɾ���ÿһ�д��һ����ɵ㣻
%��4������prf1��Ӧ�������1��sigma-profile��Ϊ2ά����column1��Ӧ����꣬column2Ϊ������profile��column3Ϊ���ֵ�profile����λΪ���λ��A^2
%��5������prf2��Ӧ�������2��sigma-profile��Ϊ2ά����column1��Ӧ����꣬column2Ϊ������profile��column3Ϊ���ֵ�profile����λΪ���λ��A^2
%��6������mixprf��Ӧ�������1�����2��ϵ�sigma-profile��Ϊ2ά����column1��Ӧ����꣬column2Ϊ������profile��column3Ϊ���ֵ�profile����λΪ���λ��A^2
%��7������Ϊ��ϵ���¶ȣ���λΪK
%��8������Ϊѡ����㽻�������ܵĲ���Ŀ���
global Aeff;
%constant for calculation
%Aeff=7.25;              %��Ч��׼Ƭ�������Ƭ��ƽ������ֵ����һ��,from the paper of COSMO-SAC(2007)
rgas=0.001987;          %���峣��kcal/mol/K
tol=0.000000001;            %������׼
numComp=length(Vcosmo); %�����
numPoint=length(X(:,1));%��ݵ���
sumA=X*(Acosmo');       %��Ϻ������,����������Ӧÿ����Ļ�����
sigma=prf1(:,1);        %sigma-profile��������
prfLevel=2;             %profile��ֳɵ�ά�ȣ�����ȡ2��ʾprofile��ֳɷ��������֣�N O F�Լ�����������H��
numSig=prfLevel*length(sigma);   %��profile�������ݵ���
%
%
%���㴿���Ƭ�εĻ��ϵ��
%     if i==1
    convPure=ones(numComp,numSig);
    seggammaPureOld=ones(numComp,numSig);
    %+���촿��ֵ�sigma-profile����
    purePrf=zeros(numSig,numComp);  %��λ��ȻΪ���λ��A^2
    purePrf(:,1)=[prf1(:,2);prf1(:,3)];
    purePrf(:,2)=[prf2(:,2);prf2(:,3)];
purePrf
    seggammaPure=ones(numComp,numSig);
    for j=1:numComp        
%         seggammaPure=ones(numComp,numSig);
        istep=0;  %��¼������
        while (sum(convPure(j,:))>=tol)
            istep=istep+1;
            convPure(j,:)=zeros(numSig,1);    %�б�Ҫ����Ϊ0��
            seggammaPureOld(j,:)=seggammaPure(j,:);
            for m=1:numSig
                summation=0.0;
                    for n=1:numSig
                        if (purePrf(n,j)~=0)
                            summation=summation+(purePrf(n,j)/Acosmo(j))*seggammaPureOld(j,n)*exp(-deltaWglobal(m,n,sigma,parameterFlag,C,T)/rgas/T);  %Ҫ��鵥λ
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
    %���㴿���Ƭ�λ��ϵ�����
%
gamma=ones(numPoint,numComp);
lngamma=zeros(numPoint,numComp);
sumgamma=zeros(numPoint,numComp);
for i=1:numPoint
    %������profile
    %+������profile��С
     TempMixprf=zeros(numSig,1);
    % mix_area=area1*X(1)+area2*X(2);%��������[��������������������]
    mixprf=prf1(:,2:3)*X(i,1)+prf2(:,2:3)*X(i,2);
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
                            summation=summation+mixprf(n)/sumA(i)*seggammaold(n)*exp(-deltaWglobal(m,n,sigma,parameterFlag,C,T)/rgas/T);
                        end
                    end                
                seggamma(m)=exp(-log(summation));
                seggamma(m)=(seggamma(m)+seggammaold(m))/2.0;
            end        
        converge=(seggamma-seggammaold).^2;
    end
    %������Ƭ�λ��ϵ�����
%
    %������ϵ��ʣ����
%     gamma=ones(1,numComp);
%     lngamma=ones(1,numComp);
%     sumgamma=ones(1,numComp);
    for k=1:numComp        
            for j=1:numSig
	%	seggamma(j)
	%	seggammaPure(k,j)
                sumgamma(i,k)=sumgamma(i,k)+((purePrf(j,k)/Aeff)*(log(seggamma(j)/seggammaPure(k,j)))); %ע���������ʽ�Ƿ���ȷ
            end        
       gamma(i,k)=exp(sumgamma(i,k));  %         gamma(i,k)=exp(sumgamma(i,k)+lngamma(i,k));��������ʽ����ȷ��
       lngamma(i,k)=log(gamma(i,k));
    end
end
end
    %���ϵ��ʣ����������
    
        
       
